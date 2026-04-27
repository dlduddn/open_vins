/*
 * OpenVINS: An Open Platform for Visual-Inertial Research
 * Copyright (C) 2018-2023 Patrick Geneva
 * Copyright (C) 2018-2023 Guoquan Huang
 * Copyright (C) 2018-2023 OpenVINS Contributors
 * Copyright (C) 2018-2019 Kevin Eckenhoff
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "VioManager.h"

#include "feat/Feature.h"
#include "feat/FeatureDatabase.h"
#include "feat/FeatureInitializer.h"
#include "track/TrackAruco.h"
#include "track/TrackDescriptor.h"
#include "track/TrackKLT.h"
#include "track/TrackSIM.h"
#include "types/Landmark.h"
#include "types/LandmarkRepresentation.h"
#include "utils/opencv_lambda_body.h"
#include "utils/print.h"
#include "utils/sensor_data.h"

#include "init/InertialInitializer.h"

#include "state/Propagator.h"
#include "state/State.h"
#include "state/StateHelper.h"
#include "update/UpdaterMSCKF.h"
#include "update/UpdaterSLAM.h"
#include "update/UpdaterZeroVelocity.h"

#include <cmath>

using namespace ov_core;
using namespace ov_type;
using namespace ov_msckf;

namespace {

template <typename Derived> bool all_finite(const Eigen::MatrixBase<Derived> &value) { return value.array().isFinite().all(); }

} // namespace

VioManager::VioManager(VioManagerOptions &params_) : thread_init_running(false), thread_init_success(false) {

  // Nice startup message
  PRINT_DEBUG("=======================================\n");
  PRINT_DEBUG("OPENVINS ON-MANIFOLD EKF IS STARTING\n");
  PRINT_DEBUG("=======================================\n");

  // Nice debug
  this->params = params_;
  params.print_and_load_estimator();
  params.print_and_load_noise();
  params.print_and_load_state();
  params.print_and_load_trackers();

  // This will globally set the thread count we will use
  // -1 will reset to the system default threading (usually the num of cores)
  cv::setNumThreads(params.num_opencv_threads);
  cv::setRNGSeed(0);

  // Create the state!!
  state = std::make_shared<State>(params.state_options);

  // Set the IMU intrinsics
  state->_calib_imu_dw->set_value(params.vec_dw);
  state->_calib_imu_dw->set_fej(params.vec_dw);
  state->_calib_imu_da->set_value(params.vec_da);
  state->_calib_imu_da->set_fej(params.vec_da);
  state->_calib_imu_tg->set_value(params.vec_tg);
  state->_calib_imu_tg->set_fej(params.vec_tg);
  state->_calib_imu_GYROtoIMU->set_value(params.q_GYROtoIMU);
  state->_calib_imu_GYROtoIMU->set_fej(params.q_GYROtoIMU);
  state->_calib_imu_ACCtoIMU->set_value(params.q_ACCtoIMU);
  state->_calib_imu_ACCtoIMU->set_fej(params.q_ACCtoIMU);

  // Timeoffset from camera to IMU
  Eigen::VectorXd temp_camimu_dt;
  temp_camimu_dt.resize(1);
  temp_camimu_dt(0) = params.calib_camimu_dt;
  state->_calib_dt_CAMtoIMU->set_value(temp_camimu_dt);
  state->_calib_dt_CAMtoIMU->set_fej(temp_camimu_dt);

  // Loop through and load each of the cameras
  state->_cam_intrinsics_cameras = params.camera_intrinsics;
  for (int i = 0; i < state->_options.num_cameras; i++) {
    state->_cam_intrinsics.at(i)->set_value(params.camera_intrinsics.at(i)->get_value());
    state->_cam_intrinsics.at(i)->set_fej(params.camera_intrinsics.at(i)->get_value());
    state->_calib_IMUtoCAM.at(i)->set_value(params.camera_extrinsics.at(i));
    state->_calib_IMUtoCAM.at(i)->set_fej(params.camera_extrinsics.at(i));
  }

  //===================================================================================
  //===================================================================================
  //===================================================================================

  // If we are recording statistics, then open our file
  if (params.record_timing_information) {
    // If the file exists, then delete it
    if (boost::filesystem::exists(params.record_timing_filepath)) {
      boost::filesystem::remove(params.record_timing_filepath);
      PRINT_INFO(YELLOW "[STATS]: found old file found, deleted...\n" RESET);
    }
    // Create the directory that we will open the file in
    boost::filesystem::path p(params.record_timing_filepath);
    boost::filesystem::create_directories(p.parent_path());
    // Open our statistics file!
    of_statistics.open(params.record_timing_filepath, std::ofstream::out | std::ofstream::app);
    // Write the header information into it
    of_statistics << "# timestamp (sec),tracking,propagation,msckf update,";
    if (state->_options.max_slam_features > 0) {
      of_statistics << "slam update,slam delayed,";
    }
    of_statistics << "re-tri & marg,total" << std::endl;
  }

  //===================================================================================
  //===================================================================================
  //===================================================================================

  // Let's make a feature extractor
  // NOTE: after we initialize we will increase the total number of feature tracks
  // NOTE: we will split the total number of features over all cameras uniformly
  int init_max_features = std::floor((double)params.init_options.init_max_features / (double)params.state_options.num_cameras);
  if (params.use_klt) {
    trackFEATS = std::shared_ptr<TrackBase>(new TrackKLT(state->_cam_intrinsics_cameras, init_max_features,
                                                         state->_options.max_aruco_features, params.use_stereo, params.histogram_method,
                                                         params.fast_threshold, params.grid_x, params.grid_y, params.min_px_dist));
  } else {
    trackFEATS = std::shared_ptr<TrackBase>(new TrackDescriptor(
        state->_cam_intrinsics_cameras, init_max_features, state->_options.max_aruco_features, params.use_stereo, params.histogram_method,
        params.fast_threshold, params.grid_x, params.grid_y, params.min_px_dist, params.knn_ratio));
  }

  // Initialize our aruco tag extractor
  if (params.use_aruco) {
    trackARUCO = std::shared_ptr<TrackBase>(new TrackAruco(state->_cam_intrinsics_cameras, state->_options.max_aruco_features,
                                                           params.use_stereo, params.histogram_method, params.downsize_aruco));
  }

  // Initialize our state propagator
  propagator = std::make_shared<Propagator>(params.imu_noises, params.gravity_mag);

  // Our state initialize
  initializer = std::make_shared<ov_init::InertialInitializer>(params.init_options, trackFEATS->get_feature_database());

  // Make the updater!
  updaterMSCKF = std::make_shared<UpdaterMSCKF>(params.msckf_options, params.featinit_options);
  updaterSLAM = std::make_shared<UpdaterSLAM>(params.slam_options, params.aruco_options, params.featinit_options);

  // If we are using zero velocity updates, then create the updater
  if (params.try_zupt) {
    updaterZUPT = std::make_shared<UpdaterZeroVelocity>(params.zupt_options, params.imu_noises, trackFEATS->get_feature_database(),
                                                        propagator, params.gravity_mag, params.zupt_max_velocity,
                                                        params.zupt_noise_multiplier, params.zupt_max_disparity);
  }
}

void VioManager::feed_measurement_imu(const ov_core::ImuData &message) {

  // The oldest time we need IMU with is the last clone
  // We shouldn't really need the whole window, but if we go backwards in time we will
  double oldest_time = state->margtimestep();
  if (oldest_time > state->_timestamp) {
    oldest_time = -1;
  }
  if (!is_initialized_vio) {
    oldest_time = message.timestamp - params.init_options.init_window_time + state->_calib_dt_CAMtoIMU->value()(0) - 0.10;
  }
  propagator->feed_imu(message, oldest_time);

  // Push back to our initializer
  if (!is_initialized_vio) {
    initializer->feed_imu(message, oldest_time);
  }

  // Push back to the zero velocity updater if it is enabled
  // No need to push back if we are just doing the zv-update at the begining and we have moved
  if (is_initialized_vio && updaterZUPT != nullptr && (!params.zupt_only_at_beginning || !has_moved_since_zupt)) {
    updaterZUPT->feed_imu(message, oldest_time);
  }
}

void VioManager::feed_measurement_simulation(double timestamp, const std::vector<int> &camids,
                                             const std::vector<std::vector<std::pair<size_t, Eigen::VectorXf>>> &feats) {

  // Start timing
  rT1 = boost::posix_time::microsec_clock::local_time();

  // Check if we actually have a simulated tracker
  // If not, recreate and re-cast the tracker to our simulation tracker
  std::shared_ptr<TrackSIM> trackSIM = std::dynamic_pointer_cast<TrackSIM>(trackFEATS);
  if (trackSIM == nullptr) {
    // Replace with the simulated tracker
    trackSIM = std::make_shared<TrackSIM>(state->_cam_intrinsics_cameras, state->_options.max_aruco_features);
    trackFEATS = trackSIM;
    // Need to also replace it in init and zv-upt since it points to the trackFEATS db pointer
    initializer = std::make_shared<ov_init::InertialInitializer>(params.init_options, trackFEATS->get_feature_database());
    if (params.try_zupt) {
      updaterZUPT = std::make_shared<UpdaterZeroVelocity>(params.zupt_options, params.imu_noises, trackFEATS->get_feature_database(),
                                                          propagator, params.gravity_mag, params.zupt_max_velocity,
                                                          params.zupt_noise_multiplier, params.zupt_max_disparity);
    }
    PRINT_WARNING(RED "[SIM]: casting our tracker to a TrackSIM object!\n" RESET);
  }

  // Feed our simulation tracker
  trackSIM->feed_measurement_simulation(timestamp, camids, feats);
  rT2 = boost::posix_time::microsec_clock::local_time();

  // Check if we should do zero-velocity, if so update the state with it
  // Note that in the case that we only use in the beginning initialization phase
  // If we have since moved, then we should never try to do a zero velocity update!
  if (is_initialized_vio && updaterZUPT != nullptr && (!params.zupt_only_at_beginning || !has_moved_since_zupt)) {
    // If the same state time, use the previous timestep decision
    if (state->_timestamp != timestamp) {
      did_zupt_update = updaterZUPT->try_update(state, timestamp);
    }
    if (did_zupt_update) {
      assert(state->_timestamp == timestamp);
      propagator->clean_old_imu_measurements(timestamp + state->_calib_dt_CAMtoIMU->value()(0) - 0.10);
      updaterZUPT->clean_old_imu_measurements(timestamp + state->_calib_dt_CAMtoIMU->value()(0) - 0.10);
      propagator->invalidate_cache();
      return;
    }
  }

  // If we do not have VIO initialization, then return an error
  if (!is_initialized_vio) {
    PRINT_ERROR(RED "[SIM]: your vio system should already be initialized before simulating features!!!\n" RESET);
    PRINT_ERROR(RED "[SIM]: initialize your system first before calling feed_measurement_simulation()!!!!\n" RESET);
    std::exit(EXIT_FAILURE);
  }

  // Call on our propagate and update function
  // Simulation is either all sync, or single camera...
  ov_core::CameraData message;
  message.timestamp = timestamp;
  for (auto const &camid : camids) {
    int width = state->_cam_intrinsics_cameras.at(camid)->w();
    int height = state->_cam_intrinsics_cameras.at(camid)->h();
    message.sensor_ids.push_back(camid);
    message.images.push_back(cv::Mat::zeros(cv::Size(width, height), CV_8UC1));
    message.masks.push_back(cv::Mat::zeros(cv::Size(width, height), CV_8UC1));
  }
  do_feature_propagate_update(message);
}

void VioManager::track_image_and_update(const ov_core::CameraData &message_const) {

  // 이 함수는 새 카메라 측정이 들어왔을 때 호출되는 메인 처리 루틴이다.
  // 큰 흐름은 image downsample -> feature tracking -> optional ZUPT ->
  // optional VIO initialization -> propagation/update 순서이다.

  // 처리 시간 측정을 시작한다. rT1/rT2는 tracking 시간 로그에 사용된다.
  rT1 = boost::posix_time::microsec_clock::local_time();

  // 입력 CameraData가 기본 조건을 만족하는지 확인한다.
  // sensor_ids가 비어 있으면 어떤 카메라에서 온 이미지인지 알 수 없다.
  assert(!message_const.sensor_ids.empty());

  // 각 camera id마다 대응되는 image가 하나씩 있어야 한다.
  assert(message_const.sensor_ids.size() == message_const.images.size());

  // 같은 CameraData 안에서 인접한 camera id가 중복되면 잘못 구성된 측정이다.
  // stereo/다중 카메라 입력에서 동일 카메라를 두 번 넣는 상황을 방지한다.
  for (size_t i = 0; i < message_const.sensor_ids.size() - 1; i++) {
    assert(message_const.sensor_ids.at(i) != message_const.sensor_ids.at(i + 1));
  }

  // 원본 입력은 const reference이므로, 필요하면 복사본을 만들어 전처리한다.
  // params.downsample_cameras가 켜져 있으면 image와 mask를 절반 해상도로 줄인다.
  ov_core::CameraData message = message_const;
  for (size_t i = 0; i < message.sensor_ids.size() && params.downsample_cameras; i++) {
    cv::Mat img = message.images.at(i);
    cv::Mat mask = message.masks.at(i);
    cv::Mat img_temp, mask_temp;

    // cv::pyrDown은 Gaussian pyramid 방식으로 해상도를 1/2로 낮춘다.
    // feature tracking 계산량을 줄이기 위한 전처리이다.
    cv::pyrDown(img, img_temp, cv::Size(img.cols / 2.0, img.rows / 2.0));
    message.images.at(i) = img_temp;

    // mask도 image와 같은 크기로 맞춰야 하므로 동일하게 downsample한다.
    cv::pyrDown(mask, mask_temp, cv::Size(mask.cols / 2.0, mask.rows / 2.0));
    message.masks.at(i) = mask_temp;
  }

  // 일반 feature tracker에 새 camera frame을 넣는다.
  // 여기서 이미지 간 feature 추적, 새 feature 검출, track database 갱신 등이 수행된다.
  trackFEATS->feed_new_camera(message);

  // VIO가 이미 초기화되었고 ARUCO tracker가 설정되어 있으면 ARUCO marker도 함께 추적한다.
  // ARUCO는 marker id를 직접 알 수 있으므로 일반 feature처럼 binocular matching을 할 필요가 적다.
  // binocular 설정에서도 내부적으로는 stereo tracking 경로를 호출한다.
  if (is_initialized_vio && trackARUCO != nullptr) {
    trackARUCO->feed_new_camera(message);
  }

  // tracking 단계가 끝난 시각을 저장한다.
  rT2 = boost::posix_time::microsec_clock::local_time();

  // Zero-velocity update(ZUPT)를 수행할 수 있는 상황인지 확인한다.
  // ZUPT는 시스템이 정지해 있다고 판단될 때 속도 0 제약으로 state를 보정하는 업데이트이다.
  // params.zupt_only_at_beginning이 켜져 있으면, 초기 정지 구간 이후 움직인 적이 있을 때는 더 이상 ZUPT를 하지 않는다.
  if (is_initialized_vio && updaterZUPT != nullptr && (!params.zupt_only_at_beginning || !has_moved_since_zupt)) {
    // 같은 timestamp에서 중복 호출된 경우에는 이전 ZUPT 판단을 재사용한다.
    // timestamp가 바뀐 경우에만 현재 시각에서 ZUPT 가능 여부를 새로 평가한다.
    if (state->_timestamp != message.timestamp) {
      did_zupt_update = updaterZUPT->try_update(state, message.timestamp);
    }

    // ZUPT가 실제로 적용되었다면, 이 frame에서는 일반 feature update를 하지 않고 종료한다.
    // ZUPT 업데이트가 state timestamp를 message timestamp까지 진행시켰는지 확인한다.
    if (did_zupt_update) {
      assert(state->_timestamp == message.timestamp);

      // 이미 사용했거나 충분히 오래된 IMU 측정은 propagation/ZUPT 버퍼에서 제거한다.
      // camera timestamp를 IMU clock으로 바꾼 뒤 약간의 여유 0.10초를 남긴다.
      propagator->clean_old_imu_measurements(message.timestamp + state->_calib_dt_CAMtoIMU->value()(0) - 0.10);
      updaterZUPT->clean_old_imu_measurements(message.timestamp + state->_calib_dt_CAMtoIMU->value()(0) - 0.10);

      // ZUPT로 state가 바뀌었으므로 propagation cache를 무효화한다.
      propagator->invalidate_cache();
      return;
    }
  }

  // 아직 VIO가 초기화되지 않았다면 현재 camera frame을 이용해 초기화를 시도한다.
  // 초기화가 성공해야 이후의 propagation/update 루틴으로 들어갈 수 있다.
  // TODO: 시스템 reset이 필요한 경우도 이 근처에서 처리할 수 있다.
  if (!is_initialized_vio) {
    is_initialized_vio = try_to_initialize(message);
    if (!is_initialized_vio) {
      // 초기화에 실패하면 이번 frame에서는 tracking 시간만 출력하고 종료한다.
      // 충분한 motion/feature/IMU 조건이 쌓일 때까지 다음 camera frame을 기다린다.
      double time_track = (rT2 - rT1).total_microseconds() * 1e-6;
      PRINT_DEBUG(BLUE "[TIME]: %.4f seconds for tracking\n" RESET, time_track);
      return;
    }
  }

  // 초기화가 완료된 상태라면, camera timestamp까지 state를 propagate하고
  // 추적된 feature/ARUCO 정보를 이용해 MSCKF, SLAM, optional update를 수행한다.
  do_feature_propagate_update(message);
}

bool VioManager::try_apply_matlab_constraint_update(const ov_core::CameraData &message) {

  // ROS wrapper(ROS1Visualizer)가 MATLAB callback을 등록하지 않았다면
  // 외부 constraint 없이 기본 OpenVINS 흐름으로만 동작한다.
  if (!matlab_constraint_callback || state == nullptr)
    return false;

  // MATLAB 외부 측정은 현재 camera frame과 같은 timestamp에 정렬되어 있다고 가정한다.
  // 따라서 직접 update할 state 변수는 camera timestamp t_k에 해당하는 IMU clone pose이다.
  auto clone_it = state->_clones_IMU.find(message.timestamp);
  if (clone_it == state->_clones_IMU.end()) {
    PRINT_WARNING(YELLOW "[MATLAB]: no IMU clone for camera timestamp %.9f; skipping external constraint\n" RESET, message.timestamp);
    return false;
  }
  const std::shared_ptr<PoseJPL> clone = clone_it->second;

  // MATLAB에 보낼 covariance와 MATLAB에서 받을 H의 column 순서를 명확히 고정한다.
  // 순서는 이 clone의 [position_error(3), orientation_error(3)]이다.
  std::vector<std::shared_ptr<Type>> pose_order = {clone->p(), clone->q()};

  // MATLAB service request에 넘길 선형화 기준 snapshot을 구성한다.
  // position/quaternion은 t_k clone의 pose이고, covariance도 같은 pose_order 기준 6x6 marginal covariance이다.
  MatlabConstraintSnapshot snapshot;
  snapshot.timestamp_cam = message.timestamp;
  snapshot.timestamp_imu = message.timestamp + state->_calib_dt_CAMtoIMU->value()(0);
  snapshot.position = clone->pos();
  snapshot.quaternion = clone->quat();
  snapshot.pose_covariance = StateHelper::get_marginal_covariance(state, pose_order);

  // NaN/Inf가 포함된 snapshot을 MATLAB으로 보내면 MATLAB 계산 또는 EKF update가 깨질 수 있다.
  // service 호출 전에 모든 숫자가 유한한지 확인한다.
  const bool snapshot_finite = std::isfinite(snapshot.timestamp_cam) && std::isfinite(snapshot.timestamp_imu) &&
                               all_finite(snapshot.position) && all_finite(snapshot.quaternion) && all_finite(snapshot.pose_covariance);
  if (!snapshot_finite) {
    PRINT_WARNING(YELLOW "[MATLAB]: non-finite clone snapshot at %.9f; skipping external constraint\n" RESET, message.timestamp);
    return false;
  }

  // 등록된 callback을 통해 ROS1Visualizer -> MATLAB service로 snapshot을 보내고,
  // MATLAB이 계산한 residual r, Jacobian H, measurement covariance R을 update에 채워 받는다.
  MatlabConstraintUpdate update;
  const bool callback_success = matlab_constraint_callback(snapshot, update);

  // callback_success=false는 service 호출 자체를 수행하지 못했거나 통신이 실패한 경우이다.
  if (!callback_success)
    return false;

  // callback은 성공했지만 MATLAB이 이 constraint를 쓰지 않겠다고 판단하면 accepted=false가 된다.
  if (!update.accepted)
    return false;

  // EKFUpdate가 요구하는 dimension을 확인한다.
  // r은 m x 1, H는 m x 6, R은 m x m이어야 한다.
  // 여기서 6은 pose_order의 [position error 3 + orientation error 3]에 대응한다.
  const int rows = static_cast<int>(update.r.rows());
  const bool dimensions_valid =
      rows > 0 && update.r.cols() == 1 && update.H.rows() == rows && update.H.cols() == 6 && update.R.rows() == rows && update.R.cols() == rows;
  if (!dimensions_valid) {
    PRINT_WARNING(YELLOW "[MATLAB]: invalid constraint dimensions (r=%dx%d, H=%dx%d, R=%dx%d); skipping update\n" RESET,
                  (int)update.r.rows(), (int)update.r.cols(), (int)update.H.rows(), (int)update.H.cols(), (int)update.R.rows(),
                  (int)update.R.cols());
    return false;
  }

  // MATLAB이 반환한 r/H/R에 NaN/Inf가 있으면 EKF update가 state와 covariance를 망가뜨릴 수 있다.
  const bool update_finite = all_finite(update.r) && all_finite(update.H) && all_finite(update.R);
  if (!update_finite) {
    PRINT_WARNING(YELLOW "[MATLAB]: non-finite H/r/R; skipping external constraint\n" RESET);
    return false;
  }

  // MATLAB constraint를 EKF update로 적용한다.
  // pose_order가 H의 column 순서와 반드시 일치해야 올바른 state block이 update된다.
  StateHelper::EKFUpdate(state, pose_order, update.H, update.r, update.R);
  PRINT_DEBUG(BLUE "[MATLAB]: applied external constraint at %.9f (%d residuals)\n" RESET, message.timestamp, rows);
  return true;
}

void VioManager::do_feature_propagate_update(const ov_core::CameraData &message) {

  //===================================================================================
  // State propagation, and clone augmentation
  //===================================================================================

  // 카메라 측정이 현재 state 시간보다 과거이면 시간 순서가 깨진 입력이다.
  // EKF propagation은 시간을 되돌릴 수 없으므로 이 frame은 사용하지 않고 종료한다.
  if (state->_timestamp > message.timestamp) {
    PRINT_WARNING(YELLOW "image received out of order, unable to do anything (prop dt = %3f)\n" RESET,
                  (message.timestamp - state->_timestamp));
    return;
  }

  // 현재 state를 camera timestamp까지 IMU propagation으로 전파하고,
  // MSCKF update에 사용할 IMU pose clone을 새로 추가한다.
  // simulation처럼 state 시간이 이미 message.timestamp와 같으면 전파할 필요가 없다.
  if (state->_timestamp != message.timestamp) {
    propagator->propagate_and_clone(state, message.timestamp);
  }
  rT3 = boost::posix_time::microsec_clock::local_time();

  // feature triangulation과 MSCKF update를 하려면 여러 시점의 clone pose가 필요하다.
  // 최소 5개 또는 max_clone_size 중 작은 값만큼 clone이 쌓일 때까지는 update하지 않고 기다린다.
  if ((int)state->_clones_IMU.size() < std::min(state->_options.max_clone_size, 5)) {
    PRINT_DEBUG("waiting for enough clone states (%d of %d)....\n", (int)state->_clones_IMU.size(),
                std::min(state->_options.max_clone_size, 5));
    // There is no visual update yet, so the best available linearization point is
    // the freshly propagated t_k clone.
    if (try_apply_matlab_constraint_update(message)) {
      propagator->invalidate_cache();
    }
    return;
  }

  // propagate_and_clone 이후에도 state 시간이 camera timestamp에 도달하지 못했다면
  // 사용할 수 있는 IMU 데이터가 부족했거나 propagation이 실패한 상황이다.
  if (state->_timestamp != message.timestamp) {
    PRINT_WARNING(RED "[PROP]: Propagator unable to propagate the state forward in time!\n" RESET);
    PRINT_WARNING(RED "[PROP]: It has been %.3f since last time we propagated\n" RESET, message.timestamp - state->_timestamp);
    return;
  }

  // 여기까지 왔다는 것은 정상적인 camera update 단계에 들어왔다는 뜻이다.
  // 이후에는 초기 정지 구간 전용 ZUPT를 더 이상 적용하지 않도록 이동 플래그를 세운다.
  has_moved_since_zupt = true;

  //===================================================================================
  // MSCKF features and KLT tracks that are SLAM features
  //===================================================================================

  // 최신 frame에서 더 이상 보이지 않는 feature들을 가져온다.
  // lost feature는 track이 끝났으므로 MSCKF update에 사용하기 좋은 후보가 된다.
  // 이미 다른 update에서 사용되어 삭제 표시된 feature는 제외한다.
  std::vector<std::shared_ptr<Feature>> feats_lost, feats_marg, feats_slam;
  feats_lost = trackFEATS->get_feature_database()->features_not_containing_newer(state->_timestamp, false, true);

  // clone window가 충분히 찼을 때는 가장 오래된 clone을 곧 marginalize해야 한다.
  // 그 clone 시각을 포함하는 feature들은 지금 update에 사용하거나 SLAM feature 후보로 넘겨야 한다.
  if ((int)state->_clones_IMU.size() > state->_options.max_clone_size || (int)state->_clones_IMU.size() > 5) {
    feats_marg = trackFEATS->get_feature_database()->features_containing(state->margtimestep(), false, true);

    // ARUCO feature는 SLAM landmark로 사용할 수 있으므로, 지연 시간 이후에만 후보로 가져온다.
    // dt_slam_delay는 초기 구간의 불안정한 landmark 추가를 막기 위한 대기 시간이다.
    if (trackARUCO != nullptr && message.timestamp - startup_time >= params.dt_slam_delay) {
      feats_slam = trackARUCO->get_feature_database()->features_containing(state->margtimestep(), false, true);
    }
  }

  // lost feature 중 현재 camera message와 관련 없는 camera stream에서만 관측된 feature는 제거한다.
  // 예를 들어 cam1 frame을 처리 중인데 cam0의 최신 frame이 아직 처리되지 않았다면,
  // cam0에서 본 feature를 너무 일찍 lost로 판단해 update에 쓰면 안 된다.
  auto it1 = feats_lost.begin();
  while (it1 != feats_lost.end()) {
    bool found_current_message_camid = false;
    for (const auto &camuvpair : (*it1)->uvs) {
      if (std::find(message.sensor_ids.begin(), message.sensor_ids.end(), camuvpair.first) != message.sensor_ids.end()) {
        found_current_message_camid = true;
        break;
      }
    }
    if (found_current_message_camid) {
      it1++;
    } else {
      it1 = feats_lost.erase(it1);
    }
  }

  // feats_lost와 feats_marg에 같은 feature가 동시에 들어가는 중복을 제거한다.
  // 이전 frame에서 lost 되었지만 marginalization timestep에도 관측이 있으면 둘 다에 잡힐 수 있다.
  it1 = feats_lost.begin();
  while (it1 != feats_lost.end()) {
    if (std::find(feats_marg.begin(), feats_marg.end(), (*it1)) != feats_marg.end()) {
      // PRINT_WARNING(YELLOW "FOUND FEATURE THAT WAS IN BOTH feats_lost and feats_marg!!!!!!\n" RESET);
      it1 = feats_lost.erase(it1);
    } else {
      it1++;
    }
  }

  // clone window 길이를 넘어선 긴 track을 찾는다.
  // 이런 feature는 MSCKF update로 소모하기 전에 SLAM landmark로 승격할 수 있는 후보가 된다.
  std::vector<std::shared_ptr<Feature>> feats_maxtracks;
  auto it2 = feats_marg.begin();
  while (it2 != feats_marg.end()) {
    // 여러 camera 중 하나라도 max_clone_size보다 긴 track을 가지고 있으면 max track으로 본다.
    bool reached_max = false;
    for (const auto &cams : (*it2)->timestamps) {
      if ((int)cams.second.size() > state->_options.max_clone_size) {
        reached_max = true;
        break;
      }
    }
    // max track이면 SLAM 후보 목록으로 옮기고, 일반 marginalization 후보에서는 제거한다.
    if (reached_max) {
      feats_maxtracks.push_back(*it2);
      it2 = feats_marg.erase(it2);
    } else {
      it2++;
    }
  }

  // 현재 state에 들어 있는 ARUCO landmark 개수를 센다.
  // 일반 SLAM feature 제한 개수와 ARUCO landmark 개수를 함께 고려하기 위해 사용된다.
  int curr_aruco_tags = 0;
  auto it0 = state->_features_SLAM.begin();
  while (it0 != state->_features_SLAM.end()) {
    if ((int)(*it0).second->_featid <= 4 * state->_options.max_aruco_features)
      curr_aruco_tags++;
    it0++;
  }

  // SLAM feature slot에 여유가 있고, startup 이후 충분한 시간이 지났으면
  // max-track 후보 중 일부를 새 SLAM feature로 추가할 준비를 한다.
  // 초기 몇 frame의 품질 낮은 landmark가 바로 들어가는 것을 dt_slam_delay로 막는다.
  if (state->_options.max_slam_features > 0 && message.timestamp - startup_time >= params.dt_slam_delay &&
      (int)state->_features_SLAM.size() < state->_options.max_slam_features + curr_aruco_tags) {
    // 현재 state에 더 추가할 수 있는 SLAM feature 수와 실제 후보 수 중 작은 값을 고른다.
    int amount_to_add = (state->_options.max_slam_features + curr_aruco_tags) - (int)state->_features_SLAM.size();
    int valid_amount = (amount_to_add > (int)feats_maxtracks.size()) ? (int)feats_maxtracks.size() : amount_to_add;

    // 추가 가능한 후보가 있으면 feats_slam으로 옮긴다.
    // 같은 feature 정보를 MSCKF와 SLAM에서 중복 사용하지 않도록 원래 후보 목록에서는 제거한다.
    if (valid_amount > 0) {
      feats_slam.insert(feats_slam.end(), feats_maxtracks.end() - valid_amount, feats_maxtracks.end());
      feats_maxtracks.erase(feats_maxtracks.end() - valid_amount, feats_maxtracks.end());
    }
  }

  // 이미 state에 들어 있는 SLAM landmark들의 최신 track을 찾아 이번 update 후보에 넣는다.
  // 현재 camera에서 더 이상 추적되지 않는 SLAM feature는 marginalize 대상으로 표시한다.
  // update 실패가 반복된 landmark도 estimator 성능을 해칠 수 있으므로 제거 대상으로 표시한다.
  for (std::pair<const size_t, std::shared_ptr<Landmark>> &landmark : state->_features_SLAM) {
    if (trackARUCO != nullptr) {
      std::shared_ptr<Feature> feat1 = trackARUCO->get_feature_database()->get_feature(landmark.second->_featid);
      if (feat1 != nullptr)
        feats_slam.push_back(feat1);
    }
    std::shared_ptr<Feature> feat2 = trackFEATS->get_feature_database()->get_feature(landmark.second->_featid);
    if (feat2 != nullptr)
      feats_slam.push_back(feat2);
    assert(landmark.second->_unique_camera_id != -1);
    bool current_unique_cam =
        std::find(message.sensor_ids.begin(), message.sensor_ids.end(), landmark.second->_unique_camera_id) != message.sensor_ids.end();
    if (feat2 == nullptr && current_unique_cam)
      landmark.second->should_marg = true;
    if (landmark.second->update_fail_count > 1)
      landmark.second->should_marg = true;
  }

  // should_marg로 표시된 오래된 SLAM landmark를 state에서 제거한다.
  // ARUCO tag landmark는 별도로 보존되도록 처리된다.
  StateHelper::marginalize_slam(state);

  // feats_slam을 이미 state에 존재하는 landmark와 새로 초기화해야 할 landmark로 나눈다.
  // 기존 landmark는 바로 EKF update에 쓰고, 새 후보는 delayed_init에서 초기화한다.
  std::vector<std::shared_ptr<Feature>> feats_slam_DELAYED, feats_slam_UPDATE;
  for (size_t i = 0; i < feats_slam.size(); i++) {
    if (state->_features_SLAM.find(feats_slam.at(i)->featid) != state->_features_SLAM.end()) {
      feats_slam_UPDATE.push_back(feats_slam.at(i));
      // PRINT_DEBUG("[UPDATE-SLAM]: found old feature %d (%d
      // measurements)\n",(int)feats_slam.at(i)->featid,(int)feats_slam.at(i)->timestamps_left.size());
    } else {
      feats_slam_DELAYED.push_back(feats_slam.at(i));
      // PRINT_DEBUG("[UPDATE-SLAM]: new feature ready %d (%d
      // measurements)\n",(int)feats_slam.at(i)->featid,(int)feats_slam.at(i)->timestamps_left.size());
    }
  }

  // SLAM update에 쓰지 않는 feature들을 MSCKF update 후보로 합친다.
  // lost, marginalization timestep 포함 feature, max-track 잔여 feature가 여기에 들어간다.
  std::vector<std::shared_ptr<Feature>> featsup_MSCKF = feats_lost;
  featsup_MSCKF.insert(featsup_MSCKF.end(), feats_marg.begin(), feats_marg.end());
  featsup_MSCKF.insert(featsup_MSCKF.end(), feats_maxtracks.begin(), feats_maxtracks.end());

  //===================================================================================
  // Now that we have a list of features, lets do the EKF update for MSCKF and SLAM!
  //===================================================================================

  // MSCKF feature 후보를 track 길이 기준으로 정렬한다.
  // 현재는 긴 track을 더 좋은 후보로 보고 뒤쪽에 남기는 단순한 선택 전략을 사용한다.
  // TODO: FOV 내 균일 분포 같은 더 좋은 feature selection 전략을 넣을 수 있다.
  auto compare_feat = [](const std::shared_ptr<Feature> &a, const std::shared_ptr<Feature> &b) -> bool {
    size_t asize = 0;
    size_t bsize = 0;
    for (const auto &pair : a->timestamps)
      asize += pair.second.size();
    for (const auto &pair : b->timestamps)
      bsize += pair.second.size();
    return asize < bsize;
  };
  std::sort(featsup_MSCKF.begin(), featsup_MSCKF.end(), compare_feat);

  // MSCKF update에 사용할 feature 수가 너무 많으면 max_msckf_in_update개만 남긴다.
  // 정렬 이후 앞쪽을 지우므로 긴 track 후보가 우선적으로 update에 사용된다.
  if ((int)featsup_MSCKF.size() > state->_options.max_msckf_in_update)
    featsup_MSCKF.erase(featsup_MSCKF.begin(), featsup_MSCKF.end() - state->_options.max_msckf_in_update);

  // 선택된 feature들로 MSCKF EKF update를 수행한다.
  updaterMSCKF->update(state, featsup_MSCKF);

  // state가 update되었으므로 이전 propagation cache는 더 이상 유효하지 않다.
  propagator->invalidate_cache();
  rT4 = boost::posix_time::microsec_clock::local_time();

  // 기존 SLAM landmark들을 EKF update한다.
  // max_slam_in_update 단위로 나누어 순차 update하면 빠르지만,
  // 모든 landmark를 한 번에 update하는 방식보다 정확도는 약간 낮을 수 있다.
  std::vector<std::shared_ptr<Feature>> feats_slam_UPDATE_TEMP;
  while (!feats_slam_UPDATE.empty()) {
    // 이번 순차 update에 사용할 feature 묶음을 잘라낸다.
    std::vector<std::shared_ptr<Feature>> featsup_TEMP;
    featsup_TEMP.insert(featsup_TEMP.begin(), feats_slam_UPDATE.begin(),
                        feats_slam_UPDATE.begin() + std::min(state->_options.max_slam_in_update, (int)feats_slam_UPDATE.size()));
    feats_slam_UPDATE.erase(feats_slam_UPDATE.begin(),
                            feats_slam_UPDATE.begin() + std::min(state->_options.max_slam_in_update, (int)feats_slam_UPDATE.size()));

    // 기존 SLAM landmark 관측으로 state와 landmark를 update한다.
    updaterSLAM->update(state, featsup_TEMP);
    feats_slam_UPDATE_TEMP.insert(feats_slam_UPDATE_TEMP.end(), featsup_TEMP.begin(), featsup_TEMP.end());
    propagator->invalidate_cache();
  }
  feats_slam_UPDATE = feats_slam_UPDATE_TEMP;
  rT5 = boost::posix_time::microsec_clock::local_time();

  // 새 SLAM feature 후보를 delayed initialization 방식으로 state에 추가한다.
  updaterSLAM->delayed_init(state, feats_slam_DELAYED);
  rT6 = boost::posix_time::microsec_clock::local_time();

  // MATLAB extension: use the post-visual-update clone at this image timestamp as the
  // linearization point, then immediately apply the returned H/r/R before cleanup/marginalization.
  if (try_apply_matlab_constraint_update(message)) {
    propagator->invalidate_cache();
  }

  //===================================================================================
  // Update our visualization feature set, and clean up the old features
  //===================================================================================

  // visualization을 위해 현재 frame에서 추적 중인 feature들을 다시 triangulate한다.
  // base camera(id 0)에서만 수행해 multi-camera 환경에서 중복 처리를 줄인다.
  if (message.sensor_ids.at(0) == 0) {

    // 현재 active track들의 3D 위치를 다시 계산한다.
    retriangulate_active_tracks(message);

    // MSCKF visualization feature 목록은 base camera 처리 시점에만 비운다.
    // 다른 camera stream에서 추가되는 feature들을 같은 update cycle 안에서 이어 붙일 수 있게 한다.
    good_features_MSCKF.clear();
  }

  // 이번 MSCKF update에 실제 사용된 feature들의 3D 위치를 visualization용으로 저장한다.
  // 사용 완료된 feature는 database cleanup에서 제거되도록 to_delete 표시를 한다.
  for (auto const &feat : featsup_MSCKF) {
    good_features_MSCKF.push_back(feat->p_FinG);
    feat->to_delete = true;
  }

  //===================================================================================
  // Cleanup, marginalize out what we don't need any more...
  //===================================================================================

  // update에 사용되어 to_delete로 표시된 feature를 tracker database에서 제거한다.
  // 이번 update에 실패한 측정은 삭제하지 않아 이후 update에서 재사용될 수 있다.
  // 새 이미지가 들어오기 전에 cleanup해야 새 측정이 잘못 삭제되지 않는다.
  trackFEATS->get_feature_database()->cleanup();
  if (trackARUCO != nullptr) {
    trackARUCO->get_feature_database()->cleanup();
  }

  // SLAM landmark의 anchor clone이 곧 marginalize될 예정이면 다른 anchor로 바꾼다.
  updaterSLAM->change_anchors(state);

  // marginalization timestep보다 오래된 feature 측정은 tracker database에서 정리한다.
  if ((int)state->_clones_IMU.size() > state->_options.max_clone_size) {
    trackFEATS->get_feature_database()->cleanup_measurements(state->margtimestep());
    if (trackARUCO != nullptr) {
      trackARUCO->get_feature_database()->cleanup_measurements(state->margtimestep());
    }
  }

  // clone window가 max_clone_size를 넘었다면 가장 오래된 IMU clone을 state에서 marginalize한다.
  StateHelper::marginalize_old_clone(state);
  rT7 = boost::posix_time::microsec_clock::local_time();

  //===================================================================================
  // Debug info, and stats tracking
  //===================================================================================

  // 각 처리 단계별 소요 시간을 계산한다.
  double time_track = (rT2 - rT1).total_microseconds() * 1e-6;
  double time_prop = (rT3 - rT2).total_microseconds() * 1e-6;
  double time_msckf = (rT4 - rT3).total_microseconds() * 1e-6;
  double time_slam_update = (rT5 - rT4).total_microseconds() * 1e-6;
  double time_slam_delay = (rT6 - rT5).total_microseconds() * 1e-6;
  double time_marg = (rT7 - rT6).total_microseconds() * 1e-6;
  double time_total = (rT7 - rT1).total_microseconds() * 1e-6;

  // debug log로 tracking, propagation, MSCKF, SLAM, marginalization 시간을 출력한다.
  PRINT_DEBUG(BLUE "[TIME]: %.4f seconds for tracking\n" RESET, time_track);
  PRINT_DEBUG(BLUE "[TIME]: %.4f seconds for propagation\n" RESET, time_prop);
  PRINT_DEBUG(BLUE "[TIME]: %.4f seconds for MSCKF update (%d feats)\n" RESET, time_msckf, (int)featsup_MSCKF.size());
  if (state->_options.max_slam_features > 0) {
    PRINT_DEBUG(BLUE "[TIME]: %.4f seconds for SLAM update (%d feats)\n" RESET, time_slam_update, (int)state->_features_SLAM.size());
    PRINT_DEBUG(BLUE "[TIME]: %.4f seconds for SLAM delayed init (%d feats)\n" RESET, time_slam_delay, (int)feats_slam_DELAYED.size());
  }
  PRINT_DEBUG(BLUE "[TIME]: %.4f seconds for re-tri & marg (%d clones in state)\n" RESET, time_marg, (int)state->_clones_IMU.size());

  std::stringstream ss;
  ss << "[TIME]: " << std::setprecision(4) << time_total << " seconds for total (camera";
  for (const auto &id : message.sensor_ids) {
    ss << " " << id;
  }
  ss << ")" << std::endl;
  PRINT_DEBUG(BLUE "%s" RESET, ss.str().c_str());

  // timing 기록 옵션이 켜져 있으면 CSV 형태로 statistics file에 저장한다.
  if (params.record_timing_information && of_statistics.is_open()) {
    // state timestamp는 camera time 기준이므로, 통계 기록에는 IMU clock 기준 timestamp를 사용한다.
    double t_ItoC = state->_calib_dt_CAMtoIMU->value()(0);
    double timestamp_inI = state->_timestamp + t_ItoC;

    // timestamp와 각 단계별 소요 시간을 한 줄로 append한다.
    of_statistics << std::fixed << std::setprecision(15) << timestamp_inI << "," << std::fixed << std::setprecision(5) << time_track << ","
                  << time_prop << "," << time_msckf << ",";
    if (state->_options.max_slam_features > 0) {
      of_statistics << time_slam_update << "," << time_slam_delay << ",";
    }
    of_statistics << time_marg << "," << time_total << std::endl;
    of_statistics.flush();
  }

  // 이전 update clone과 현재 IMU 위치 차이를 누적해 이동 거리를 추정한다.
  if (timelastupdate != -1 && state->_clones_IMU.find(timelastupdate) != state->_clones_IMU.end()) {
    Eigen::Matrix<double, 3, 1> dx = state->_imu->pos() - state->_clones_IMU.at(timelastupdate)->pos();
    distance += dx.norm();
  }
  timelastupdate = message.timestamp;

  // 현재 IMU pose와 bias를 콘솔에 출력한다.
  PRINT_INFO("q_GtoI = %.3f,%.3f,%.3f,%.3f | p_IinG = %.3f,%.3f,%.3f | dist = %.2f (meters)\n", state->_imu->quat()(0),
             state->_imu->quat()(1), state->_imu->quat()(2), state->_imu->quat()(3), state->_imu->pos()(0), state->_imu->pos()(1),
             state->_imu->pos()(2), distance);
  PRINT_INFO("bg = %.4f,%.4f,%.4f | ba = %.4f,%.4f,%.4f\n", state->_imu->bias_g()(0), state->_imu->bias_g()(1), state->_imu->bias_g()(2),
             state->_imu->bias_a()(0), state->_imu->bias_a()(1), state->_imu->bias_a()(2));

  // camera-IMU time offset을 calibration 중이면 현재 추정값을 출력한다.
  if (state->_options.do_calib_camera_timeoffset) {
    PRINT_INFO("camera-imu timeoffset = %.5f\n", state->_calib_dt_CAMtoIMU->value()(0));
  }

  // camera intrinsics를 calibration 중이면 각 카메라의 현재 추정값을 출력한다.
  if (state->_options.do_calib_camera_intrinsics) {
    for (int i = 0; i < state->_options.num_cameras; i++) {
      std::shared_ptr<Vec> calib = state->_cam_intrinsics.at(i);
      PRINT_INFO("cam%d intrinsics = %.3f,%.3f,%.3f,%.3f | %.3f,%.3f,%.3f,%.3f\n", (int)i, calib->value()(0), calib->value()(1),
                 calib->value()(2), calib->value()(3), calib->value()(4), calib->value()(5), calib->value()(6), calib->value()(7));
    }
  }

  // IMU-to-camera extrinsics를 calibration 중이면 각 카메라의 pose 추정값을 출력한다.
  if (state->_options.do_calib_camera_pose) {
    for (int i = 0; i < state->_options.num_cameras; i++) {
      std::shared_ptr<PoseJPL> calib = state->_calib_IMUtoCAM.at(i);
      PRINT_INFO("cam%d extrinsics = %.3f,%.3f,%.3f,%.3f | %.3f,%.3f,%.3f\n", (int)i, calib->quat()(0), calib->quat()(1), calib->quat()(2),
                 calib->quat()(3), calib->pos()(0), calib->pos()(1), calib->pos()(2));
    }
  }

  // IMU intrinsics calibration 값을 출력한다.
  // Kalibr 모델과 RPNG 모델은 파라미터 구조가 달라 출력 형식이 다르다.
  if (state->_options.do_calib_imu_intrinsics && state->_options.imu_model == StateOptions::ImuModel::KALIBR) {
    PRINT_INFO("q_GYROtoI = %.3f,%.3f,%.3f,%.3f\n", state->_calib_imu_GYROtoIMU->value()(0), state->_calib_imu_GYROtoIMU->value()(1),
               state->_calib_imu_GYROtoIMU->value()(2), state->_calib_imu_GYROtoIMU->value()(3));
  }
  if (state->_options.do_calib_imu_intrinsics && state->_options.imu_model == StateOptions::ImuModel::RPNG) {
    PRINT_INFO("q_ACCtoI = %.3f,%.3f,%.3f,%.3f\n", state->_calib_imu_ACCtoIMU->value()(0), state->_calib_imu_ACCtoIMU->value()(1),
               state->_calib_imu_ACCtoIMU->value()(2), state->_calib_imu_ACCtoIMU->value()(3));
  }
  if (state->_options.do_calib_imu_intrinsics && state->_options.imu_model == StateOptions::ImuModel::KALIBR) {
    PRINT_INFO("Dw = | %.4f,%.4f,%.4f | %.4f,%.4f | %.4f |\n", state->_calib_imu_dw->value()(0), state->_calib_imu_dw->value()(1),
               state->_calib_imu_dw->value()(2), state->_calib_imu_dw->value()(3), state->_calib_imu_dw->value()(4),
               state->_calib_imu_dw->value()(5));
    PRINT_INFO("Da = | %.4f,%.4f,%.4f | %.4f,%.4f | %.4f |\n", state->_calib_imu_da->value()(0), state->_calib_imu_da->value()(1),
               state->_calib_imu_da->value()(2), state->_calib_imu_da->value()(3), state->_calib_imu_da->value()(4),
               state->_calib_imu_da->value()(5));
  }
  if (state->_options.do_calib_imu_intrinsics && state->_options.imu_model == StateOptions::ImuModel::RPNG) {
    PRINT_INFO("Dw = | %.4f | %.4f,%.4f | %.4f,%.4f,%.4f |\n", state->_calib_imu_dw->value()(0), state->_calib_imu_dw->value()(1),
               state->_calib_imu_dw->value()(2), state->_calib_imu_dw->value()(3), state->_calib_imu_dw->value()(4),
               state->_calib_imu_dw->value()(5));
    PRINT_INFO("Da = | %.4f | %.4f,%.4f | %.4f,%.4f,%.4f |\n", state->_calib_imu_da->value()(0), state->_calib_imu_da->value()(1),
               state->_calib_imu_da->value()(2), state->_calib_imu_da->value()(3), state->_calib_imu_da->value()(4),
               state->_calib_imu_da->value()(5));
  }
  if (state->_options.do_calib_imu_intrinsics && state->_options.do_calib_imu_g_sensitivity) {
    PRINT_INFO("Tg = | %.4f,%.4f,%.4f |  %.4f,%.4f,%.4f | %.4f,%.4f,%.4f |\n", state->_calib_imu_tg->value()(0),
               state->_calib_imu_tg->value()(1), state->_calib_imu_tg->value()(2), state->_calib_imu_tg->value()(3),
               state->_calib_imu_tg->value()(4), state->_calib_imu_tg->value()(5), state->_calib_imu_tg->value()(6),
               state->_calib_imu_tg->value()(7), state->_calib_imu_tg->value()(8));
  }
}
