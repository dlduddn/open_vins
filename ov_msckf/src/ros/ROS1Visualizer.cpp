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

#include "ROS1Visualizer.h"

#include "core/VioManager.h"
#include "ov_msckf/PoseSnapshotToMatlab.h"
#include "ros/ROSVisualizerHelper.h"
#include "sim/Simulator.h"
#include "state/Propagator.h"
#include "state/State.h"
#include "state/StateHelper.h"
#include "utils/dataset_reader.h"
#include "utils/print.h"
#include "utils/sensor_data.h"

#include <cmath>

using namespace ov_core;
using namespace ov_type;
using namespace ov_msckf;

// -----------------------------------------------------------------------------
// Local helpers
// -----------------------------------------------------------------------------
// Small file-local utilities used by the MATLAB snapshot extension.
namespace {

template <typename Derived> bool all_finite(const Eigen::MatrixBase<Derived> &value) { return value.array().isFinite().all(); }

} // namespace

// -----------------------------------------------------------------------------
// Construction: ROS publishers, optional MATLAB client, and output files
// -----------------------------------------------------------------------------
ROS1Visualizer::ROS1Visualizer(std::shared_ptr<ros::NodeHandle> nh, std::shared_ptr<VioManager> app, std::shared_ptr<Simulator> sim)
    : _nh(nh), _app(app), _sim(sim), thread_update_running(false) {

  // 생성자에서는 ROS visualization에 필요한 publisher와 optional 기능을 한 번에 준비한다.
  // nh는 ROS node handle, app은 VIO core(VioManager), sim은 simulation 사용 시 groundtruth를 제공한다.

  // MATLAB extension 설정값을 VioManagerOptions에서 읽어 ROS1Visualizer 멤버 변수에 저장한다.
  // 이후 MATLAB service client 생성과 callback 요청에서 이 값을 사용한다.
  const auto params = _app->get_params();
  matlab_snapshot_enable = params.matlab_snapshot_enable;
  matlab_snapshot_service_name = params.matlab_snapshot_service_name;
  matlab_snapshot_timeout_sec = params.matlab_snapshot_timeout_sec;
  matlab_snapshot_debug_log = params.matlab_snapshot_debug_log;

  // TF broadcaster를 만든다.
  // publish_global_to_imu_tf 또는 calibration tf 옵션이 켜져 있으면 transform을 ROS TF로 publish할 때 사용된다.
  mTfBr = std::make_shared<tf::TransformBroadcaster>();

  // image_transport는 ROS image topic을 publish하기 위한 전용 helper이다.
  // 일반 ROS publisher보다 image transport plugin/compression 경로와 잘 맞는다.
  image_transport::ImageTransport it(*_nh);

  // 아래 publisher들은 OpenVINS의 기본 출력이다.
  // state estimate, path, feature cloud, tracking image, groundtruth, loop-closure 정보를 ROS topic으로 내보낸다.

  // 현재 IMU pose/covariance, odometry, 누적 path를 publish할 topic을 만든다.
  pub_poseimu = nh->advertise<geometry_msgs::PoseWithCovarianceStamped>("poseimu", 2);
  PRINT_DEBUG("Publishing: %s\n", pub_poseimu.getTopic().c_str());
  pub_odomimu = nh->advertise<nav_msgs::Odometry>("odomimu", 2);
  PRINT_DEBUG("Publishing: %s\n", pub_odomimu.getTopic().c_str());
  pub_pathimu = nh->advertise<nav_msgs::Path>("pathimu", 2);
  PRINT_DEBUG("Publishing: %s\n", pub_pathimu.getTopic().c_str());

  // 추정/관리 중인 3D feature point들을 종류별로 publish한다.
  // MSCKF feature, SLAM landmark, ARUCO marker, simulator feature를 RViz 등에서 볼 수 있다.
  pub_points_msckf = nh->advertise<sensor_msgs::PointCloud2>("points_msckf", 2);
  PRINT_DEBUG("Publishing: %s\n", pub_points_msckf.getTopic().c_str());
  pub_points_slam = nh->advertise<sensor_msgs::PointCloud2>("points_slam", 2);
  PRINT_DEBUG("Publishing: %s\n", pub_points_msckf.getTopic().c_str());
  pub_points_aruco = nh->advertise<sensor_msgs::PointCloud2>("points_aruco", 2);
  PRINT_DEBUG("Publishing: %s\n", pub_points_aruco.getTopic().c_str());
  pub_points_sim = nh->advertise<sensor_msgs::PointCloud2>("points_sim", 2);
  PRINT_DEBUG("Publishing: %s\n", pub_points_sim.getTopic().c_str());

  // feature tracking 결과를 그린 image를 publish한다.
  // track history나 feature movement를 시각적으로 확인하는 용도이다.
  it_pub_tracks = it.advertise("trackhist", 2);
  PRINT_DEBUG("Publishing: %s\n", it_pub_tracks.getTopic().c_str());

  // groundtruth pose/path가 있으면 publish할 topic을 만든다.
  // dataset groundtruth 또는 simulator groundtruth와 estimate를 비교할 때 사용된다.
  pub_posegt = nh->advertise<geometry_msgs::PoseStamped>("posegt", 2);
  PRINT_DEBUG("Publishing: %s\n", pub_posegt.getTopic().c_str());
  pub_pathgt = nh->advertise<nav_msgs::Path>("pathgt", 2);
  PRINT_DEBUG("Publishing: %s\n", pub_pathgt.getTopic().c_str());

  // loop-closure 관련 pose, feature, calibration, depth image 출력을 위한 publisher들이다.
  // loop closure 모듈이 켜져 있을 때 후단 모듈이나 시각화 도구가 이 topic들을 구독할 수 있다.
  pub_loop_pose = nh->advertise<nav_msgs::Odometry>("loop_pose", 2);
  pub_loop_point = nh->advertise<sensor_msgs::PointCloud>("loop_feats", 2);
  pub_loop_extrinsic = nh->advertise<nav_msgs::Odometry>("loop_extrinsic", 2);
  pub_loop_intrinsics = nh->advertise<sensor_msgs::CameraInfo>("loop_intrinsics", 2);
  it_pub_loop_img_depth = it.advertise("loop_depth", 2);
  it_pub_loop_img_depth_color = it.advertise("loop_depth_colored", 2);

  // MATLAB extension이 켜져 있고 service 이름이 설정되어 있으면 MATLAB service client를 만든다.
  // 이 client는 VioManager core가 MATLAB constraint를 요청할 때 사용된다.
  if (matlab_snapshot_enable && !matlab_snapshot_service_name.empty()) {
    matlab_snapshot_client = _nh->serviceClient<ov_msckf::PoseSnapshotToMatlab>(matlab_snapshot_service_name, false);

    // VioManager는 ROS를 직접 알지 않도록 callback 인터페이스만 가진다.
    // 여기서 ROS1Visualizer가 람다를 등록해, core가 snapshot/update를 넘기면
    // request_matlab_constraint_update()를 통해 실제 ROS service call을 수행하게 한다.
    // [this]는 현재 ROS1Visualizer 객체의 멤버 함수와 service client에 접근하기 위한 캡처이다.
    _app->set_matlab_constraint_callback(
        [this](const MatlabConstraintSnapshot &snapshot, MatlabConstraintUpdate &update) { return request_matlab_constraint_update(snapshot, update); });

    // debug log가 켜져 있으면 어떤 MATLAB service에 연결하도록 설정되었는지 출력한다.
    if (matlab_snapshot_debug_log) {
      PRINT_DEBUG("Configured MATLAB snapshot service client: %s (timeout=%.3f sec)\n", matlab_snapshot_service_name.c_str(),
                  matlab_snapshot_timeout_sec);
    }
  } else if (matlab_snapshot_enable && matlab_snapshot_debug_log) {
    // 기능은 켜졌지만 service 이름이 비어 있으면 실제 요청은 보낼 수 없다.
    PRINT_DEBUG("MATLAB snapshot requests enabled but service name is empty; requests will be skipped.\n");
  }

  // ROS parameter에서 TF publish 옵션을 읽는다.
  // global->IMU transform과 calibration transform을 TF tree에 올릴지 결정한다.
  nh->param<bool>("publish_global_to_imu_tf", publish_global2imu_tf, true);
  nh->param<bool>("publish_calibration_tf", publish_calibration_tf, true);

  // simulation이 아닌 dataset 실행에서 path_gt parameter가 있으면 groundtruth CSV를 로드한다.
  // 파일은 ASL dataset 형식의 CSV여야 한다.
  if (nh->hasParam("path_gt") && _sim == nullptr) {
    std::string path_to_gt;
    nh->param<std::string>("path_gt", path_to_gt, "");
    if (!path_to_gt.empty()) {
      DatasetReader::load_gt_file(path_to_gt, gt_states);
      PRINT_DEBUG("gt file path is: %s\n", path_to_gt.c_str());
    }
  }

  // 전체 state 추정값과 표준편차를 파일로 저장할지 ROS parameter에서 읽는다.
  // 저장이 켜져 있으면 출력 경로를 준비하고 파일 header를 쓴다.
  nh->param<bool>("save_total_state", save_total_state, false);
  if (save_total_state) {

    // 저장할 estimate, standard deviation, groundtruth 파일 경로를 읽는다.
    std::string filepath_est, filepath_std, filepath_gt;
    nh->param<std::string>("filepath_est", filepath_est, "state_estimate.txt");
    nh->param<std::string>("filepath_std", filepath_std, "state_deviation.txt");
    nh->param<std::string>("filepath_gt", filepath_gt, "state_groundtruth.txt");

    // 같은 경로의 이전 결과 파일이 있으면 새 실행 결과와 섞이지 않도록 삭제한다.
    if (boost::filesystem::exists(filepath_est))
      boost::filesystem::remove(filepath_est);
    if (boost::filesystem::exists(filepath_std))
      boost::filesystem::remove(filepath_std);

    // 출력 파일이 들어갈 폴더가 없으면 생성한다.
    boost::filesystem::create_directories(boost::filesystem::path(filepath_est.c_str()).parent_path());
    boost::filesystem::create_directories(boost::filesystem::path(filepath_std.c_str()).parent_path());

    // estimate/std 파일을 열고 column 설명 header를 기록한다.
    of_state_est.open(filepath_est.c_str());
    of_state_std.open(filepath_std.c_str());
    of_state_est << "# timestamp(s) q p v bg ba cam_imu_dt num_cam cam0_k cam0_d cam0_rot cam0_trans ... imu_model dw da tg wtoI atoI etc"
                 << std::endl;
    of_state_std << "# timestamp(s) q p v bg ba cam_imu_dt num_cam cam0_k cam0_d cam0_rot cam0_trans ... imu_model dw da tg wtoI atoI etc"
                 << std::endl;

    // simulation 실행에서는 simulator groundtruth도 별도 파일로 저장한다.
    if (_sim != nullptr) {
      if (boost::filesystem::exists(filepath_gt))
        boost::filesystem::remove(filepath_gt);
      boost::filesystem::create_directories(boost::filesystem::path(filepath_gt.c_str()).parent_path());
      of_state_gt.open(filepath_gt.c_str());
      of_state_gt << "# timestamp(s) q p v bg ba cam_imu_dt num_cam cam0_k cam0_d cam0_rot cam0_trans ... imu_model dw da tg wtoI atoI etc"
                  << std::endl;
    }
  }

  // image publish를 별도 thread에서 돌릴지 설정값에 따라 결정한다.
  // multi-thread publisher가 켜져 있으면 20Hz로 publish_images()를 반복 호출한다.
  if (_app->get_params().use_multi_threading_pubs) {
    std::thread thread([&] {
      ros::Rate loop_rate(20);
      while (ros::ok()) {
        publish_images();
        loop_rate.sleep();
      }
    });
    thread.detach();
  }
}

// -----------------------------------------------------------------------------
// MATLAB extension: request a linearized constraint for a specific clone snapshot
// -----------------------------------------------------------------------------
bool ROS1Visualizer::request_matlab_constraint_update(const MatlabConstraintSnapshot &snapshot, MatlabConstraintUpdate &update) {

  // 설정에서 MATLAB constraint 기능이 꺼져 있으면 estimator 쪽에는 "호출하지 못함"으로 알려준다.
  // 이 경우 VioManager는 외부 constraint 없이 기존 OpenVINS update만 수행한다.
  if (!matlab_snapshot_enable)
    return false;

  // service 이름이 비어 있으면 ROS service client가 연결할 대상이 없다.
  // false를 반환하면 VioManager는 MATLAB update를 적용하지 않는다.
  if (matlab_snapshot_service_name.empty()) {
    if (matlab_snapshot_debug_log) {
      PRINT_DEBUG("Skipping MATLAB constraint request because service name is empty.\n");
    } else {
      ROS_WARN_STREAM_THROTTLE(5.0, "Skipping MATLAB constraint request because matlab_snapshot_service_name is empty.");
    }
    return false;
  }

  // VioManager가 넘겨준 snapshot은 MATLAB이 선형화 기준점으로 사용할 pose이다.
  // timestamp, position, quaternion, covariance 중 NaN/Inf가 있으면 MATLAB 계산과 EKF update가
  // 모두 위험해지므로 service 호출 전에 차단한다.
  const bool values_finite = std::isfinite(snapshot.timestamp_cam) && std::isfinite(snapshot.timestamp_imu) && all_finite(snapshot.position) &&
                             all_finite(snapshot.quaternion) && all_finite(snapshot.pose_covariance);
  if (!values_finite) {
    if (matlab_snapshot_debug_log) {
      PRINT_DEBUG("Skipping MATLAB constraint request due to non-finite snapshot values: t_cam=%.6f t_imu=%.6f\n", snapshot.timestamp_cam,
                  snapshot.timestamp_imu);
    } else {
      ROS_WARN_STREAM_THROTTLE(5.0, "Skipping MATLAB constraint request due to non-finite snapshot values.");
    }
    return false;
  }

  // ROS service request 메시지를 만든다.
  // 여기 담기는 pose는 현재 IMU state가 아니라, VioManager가 선택한 t_k clone pose이다.
  ov_msckf::PoseSnapshotToMatlab snapshot_srv;
  snapshot_srv.request.timestamp_cam = snapshot.timestamp_cam;
  snapshot_srv.request.timestamp_imu = snapshot.timestamp_imu;

  // Eigen vector를 ROS 고정 길이 배열 필드로 복사한다.
  // quaternion은 OpenVINS 내부 convention인 JPL [x y z w] 순서 그대로 보낸다.
  for (int i = 0; i < 3; i++) {
    snapshot_srv.request.position[i] = snapshot.position(i);
  }
  for (int i = 0; i < 4; i++) {
    snapshot_srv.request.quaternion[i] = snapshot.quaternion(i);
  }

  // covariance는 MATLAB에서 reshape하기 쉽도록 row-major 1차원 배열로 펼친다.
  // 이 서비스의 6D pose error 순서는 [position_error(3), orientation_error(3)]이다.
  // MATLAB이 반환하는 H의 6개 column도 반드시 같은 순서를 따라야 한다.
  for (int r = 0; r < 6; r++) {
    for (int c = 0; c < 6; c++) {
      snapshot_srv.request.pose_covariance_row_major[6 * r + c] = snapshot.pose_covariance(r, c);
    }
  }

  // MATLAB service server가 ROS master에 등록되어 있는지 확인한다.
  // timeout이 양수이면 해당 시간만큼 등록을 기다리고, 0 이하이면 현재 존재 여부만 즉시 확인한다.
  // 주의: 이 timeout은 service 등록 대기 시간이고, 아래 call() 자체의 실행 시간 제한은 아니다.
  const bool service_available = (matlab_snapshot_timeout_sec > 0.0)
                                     ? matlab_snapshot_client.waitForExistence(ros::Duration(matlab_snapshot_timeout_sec))
                                     : matlab_snapshot_client.exists();
  if (!service_available) {
    if (matlab_snapshot_debug_log) {
      PRINT_DEBUG("MATLAB constraint service unavailable: %s (timeout=%.3f sec)\n", matlab_snapshot_service_name.c_str(),
                  matlab_snapshot_timeout_sec);
    } else {
      ROS_WARN_STREAM_THROTTLE(5.0, "Skipping MATLAB constraint request because service '" << matlab_snapshot_service_name
                                                                                            << "' is unavailable (timeout="
                                                                                            << matlab_snapshot_timeout_sec << " sec).");
    }
    return false;
  }

  // 동기식 service 호출이다.
  // MATLAB callback이 H/r/R을 계산해서 response를 반환할 때까지 camera update thread가 여기서 대기한다.
  // 따라서 이 함수가 반환된 직후 VioManager가 같은 선형화 기준점에 EKF update를 걸 수 있다.
  const bool call_success = matlab_snapshot_client.call(snapshot_srv);

  // call_success=false는 ROS service transport 자체가 실패한 경우이다.
  // MATLAB이 계산은 했지만 constraint를 쓰지 않겠다고 판단한 경우는 accepted=false로 구분한다.
  if (!call_success) {
    if (!matlab_snapshot_debug_log) {
      ROS_WARN_STREAM_THROTTLE(5.0, "Skipping MATLAB constraint request because service call to '" << matlab_snapshot_service_name
                                                                                                     << "' failed.");
    }
    return false;
  }

  // MATLAB의 판단 결과를 VioManager에 넘긴다.
  // 이 함수의 return true는 "service call은 정상 처리됨"이고,
  // update.accepted가 실제 EKF update 적용 여부를 결정한다.
  update.accepted = snapshot_srv.response.accepted;
  update.status_message = snapshot_srv.response.status_message;
  if (!update.accepted)
    return true;

  // MATLAB이 accepted=true를 보냈다면 residual/Jacobian/noise payload를 Eigen 행렬로 복원한다.
  // residual_rows = m, H = m x 6, R = m x m 형태가 되어야 한다.
  const int residual_rows = static_cast<int>(snapshot_srv.response.residual.size());
  const int jacobian_rows = static_cast<int>(snapshot_srv.response.jacobian_rows);
  const int jacobian_cols = static_cast<int>(snapshot_srv.response.jacobian_cols);
  const int measurement_cov_rows = static_cast<int>(snapshot_srv.response.measurement_cov_rows);
  const int measurement_cov_cols = static_cast<int>(snapshot_srv.response.measurement_cov_cols);

  // response의 matrix dimension 필드와 실제 1차원 배열 길이가 일치하는지 확인한다.
  // 이 검사를 하지 않으면 아래 at()에서 예외가 나거나, 더 나쁘게는 잘못된 H/R로 EKF update가 들어갈 수 있다.
  const size_t expected_jacobian_size = static_cast<size_t>(jacobian_rows) * static_cast<size_t>(jacobian_cols);
  const size_t expected_measurement_cov_size = static_cast<size_t>(measurement_cov_rows) * static_cast<size_t>(measurement_cov_cols);

  if (snapshot_srv.response.jacobian_row_major.size() != expected_jacobian_size ||
      snapshot_srv.response.measurement_cov_row_major.size() != expected_measurement_cov_size) {
    update.accepted = false;
    if (matlab_snapshot_debug_log) {
      PRINT_DEBUG("MATLAB constraint rejected due to malformed matrix payloads: H=%u x %u has %zu values, R=%u x %u has %zu values\n",
                  snapshot_srv.response.jacobian_rows, snapshot_srv.response.jacobian_cols,
                  snapshot_srv.response.jacobian_row_major.size(),
                  snapshot_srv.response.measurement_cov_rows, snapshot_srv.response.measurement_cov_cols,
                  snapshot_srv.response.measurement_cov_row_major.size());
    }
    return true;
  }

  // residual r을 Eigen::VectorXd로 복사한다.
  // StateHelper::EKFUpdate()는 dx = K * res 형태를 쓰므로 residual 부호 convention을 MATLAB과 맞춰야 한다.
  update.r.resize(residual_rows);
  for (int r = 0; r < residual_rows; r++) {
    update.r(r) = snapshot_srv.response.residual.at(r);
  }

  // MATLAB이 row-major로 펼쳐 보낸 H를 Eigen 행렬로 복원한다.
  // H column 순서는 VioManager의 H_order와 같아야 한다:
  // [position_error(3), orientation_error(3)] for the t_k clone.
  update.H.resize(jacobian_rows, jacobian_cols);
  for (int r = 0; r < jacobian_rows; r++) {
    for (int c = 0; c < jacobian_cols; c++) {
      update.H(r, c) = snapshot_srv.response.jacobian_row_major.at(jacobian_cols * r + c);
    }
  }

  // MATLAB이 row-major로 펼쳐 보낸 measurement noise covariance R을 복원한다.
  // R은 residual dimension과 같은 정방행렬이어야 하며, 최종 dimension 검사는 VioManager에서 한 번 더 수행한다.
  update.R.resize(measurement_cov_rows, measurement_cov_cols);
  for (int r = 0; r < measurement_cov_rows; r++) {
    for (int c = 0; c < measurement_cov_cols; c++) {
      update.R(r, c) = snapshot_srv.response.measurement_cov_row_major.at(measurement_cov_cols * r + c);
    }
  }

  // 여기까지 오면 service 호출과 payload 복원이 성공했다는 뜻이다.
  // 실제 EKF update 적용 여부는 VioManager가 update.accepted와 dimension/finite 검사를 보고 결정한다.
  return true;
}

// -----------------------------------------------------------------------------
// Input setup: subscribe to IMU and camera ROS topics
// -----------------------------------------------------------------------------
void ROS1Visualizer::setup_subscribers(std::shared_ptr<ov_core::YamlParser> parser) {

  // We need a valid parser
  assert(parser != nullptr);

  // Create imu subscriber (handle legacy ros param info)
  std::string topic_imu;
  _nh->param<std::string>("topic_imu", topic_imu, "/imu0");
  parser->parse_external("relative_config_imu", "imu0", "rostopic", topic_imu);
  sub_imu = _nh->subscribe(topic_imu, 1000, &ROS1Visualizer::callback_inertial, this);
  PRINT_INFO("subscribing to IMU: %s\n", topic_imu.c_str());

  // Logic for sync stereo subscriber
  // https://answers.ros.org/question/96346/subscribe-to-two-image_raws-with-one-function/?answer=96491#post-id-96491
  if (_app->get_params().state_options.num_cameras == 2) {
    // Read in the topics
    std::string cam_topic0, cam_topic1;
    _nh->param<std::string>("topic_camera" + std::to_string(0), cam_topic0, "/cam" + std::to_string(0) + "/image_raw");
    _nh->param<std::string>("topic_camera" + std::to_string(1), cam_topic1, "/cam" + std::to_string(1) + "/image_raw");
    parser->parse_external("relative_config_imucam", "cam" + std::to_string(0), "rostopic", cam_topic0);
    parser->parse_external("relative_config_imucam", "cam" + std::to_string(1), "rostopic", cam_topic1);
    // Create sync filter (they have unique pointers internally, so we have to use move logic here...)
    auto image_sub0 = std::make_shared<message_filters::Subscriber<sensor_msgs::Image>>(*_nh, cam_topic0, 1);
    auto image_sub1 = std::make_shared<message_filters::Subscriber<sensor_msgs::Image>>(*_nh, cam_topic1, 1);
    auto sync = std::make_shared<message_filters::Synchronizer<sync_pol>>(sync_pol(10), *image_sub0, *image_sub1);
    sync->registerCallback(boost::bind(&ROS1Visualizer::callback_stereo, this, _1, _2, 0, 1));
    // Append to our vector of subscribers
    sync_cam.push_back(sync);
    sync_subs_cam.push_back(image_sub0);
    sync_subs_cam.push_back(image_sub1);
    PRINT_INFO("subscribing to cam (stereo): %s\n", cam_topic0.c_str());
    PRINT_INFO("subscribing to cam (stereo): %s\n", cam_topic1.c_str());
  } else {
    // Now we should add any non-stereo callbacks here
    for (int i = 0; i < _app->get_params().state_options.num_cameras; i++) {
      // read in the topic
      std::string cam_topic;
      _nh->param<std::string>("topic_camera" + std::to_string(i), cam_topic, "/cam" + std::to_string(i) + "/image_raw");
      parser->parse_external("relative_config_imucam", "cam" + std::to_string(i), "rostopic", cam_topic);
      // create subscriber
      subs_cam.push_back(_nh->subscribe<sensor_msgs::Image>(cam_topic, 10, boost::bind(&ROS1Visualizer::callback_monocular, this, _1, i)));
      PRINT_INFO("subscribing to cam (mono): %s\n", cam_topic.c_str());
    }
  }
}

// -----------------------------------------------------------------------------
// Output cycle: publish the current estimator state and visualization products
// -----------------------------------------------------------------------------
void ROS1Visualizer::visualize() {

  // Return if we have already visualized
  if (last_visualization_timestamp == _app->get_state()->_timestamp && _app->initialized())
    return;
  last_visualization_timestamp = _app->get_state()->_timestamp;

  // Start timing
  // boost::posix_time::ptime rT0_1, rT0_2;
  // rT0_1 = boost::posix_time::microsec_clock::local_time();

  // publish current image (only if not multi-threaded)
  if (!_app->get_params().use_multi_threading_pubs)
    publish_images();

  // Return if we have not inited
  if (!_app->initialized())
    return;

  // Save the start time of this dataset
  if (!start_time_set) {
    rT1 = boost::posix_time::microsec_clock::local_time();
    start_time_set = true;
  }

  // publish state
  publish_state();

  // publish points
  publish_features();

  // Publish gt if we have it
  publish_groundtruth();

  // Publish keyframe information
  publish_loopclosure_information();

  // Save total state
  if (save_total_state) {
    ROSVisualizerHelper::sim_save_total_state_to_file(_app->get_state(), _sim, of_state_est, of_state_std, of_state_gt);
  }

  // Print how much time it took to publish / displaying things
  // rT0_2 = boost::posix_time::microsec_clock::local_time();
  // double time_total = (rT0_2 - rT0_1).total_microseconds() * 1e-6;
  // PRINT_DEBUG(BLUE "[TIME]: %.4f seconds for visualization\n" RESET, time_total);
}

// -----------------------------------------------------------------------------
// High-rate odometry: propagate the latest state to each IMU timestamp
// -----------------------------------------------------------------------------
void ROS1Visualizer::visualize_odometry(double timestamp) {

  // Return if we have not inited
  if (!_app->initialized())
    return;

  // Get fast propagate state at the desired timestamp
  std::shared_ptr<State> state = _app->get_state();
  Eigen::Matrix<double, 13, 1> state_plus = Eigen::Matrix<double, 13, 1>::Zero();
  Eigen::Matrix<double, 12, 12> cov_plus = Eigen::Matrix<double, 12, 12>::Zero();
  if (!_app->get_propagator()->fast_state_propagate(state, timestamp, state_plus, cov_plus))
    return;

  //  // Get the simulated groundtruth so we can evaulate the error in respect to it
  //  // NOTE: we get the true time in the IMU clock frame
  //  if (_sim != nullptr) {
  //    Eigen::Matrix<double, 17, 1> state_gt;
  //    if (_sim->get_state(timestamp, state_gt)) {
  //      // Difference between positions
  //      double dx = state_plus(4, 0) - state_gt(5, 0);
  //      double dy = state_plus(5, 0) - state_gt(6, 0);
  //      double dz = state_plus(6, 0) - state_gt(7, 0);
  //      double err_pos = std::sqrt(dx * dx + dy * dy + dz * dz);
  //      // Quaternion error
  //      Eigen::Matrix<double, 4, 1> quat_gt, quat_st, quat_diff;
  //      quat_gt << state_gt(1, 0), state_gt(2, 0), state_gt(3, 0), state_gt(4, 0);
  //      quat_st << state_plus(0, 0), state_plus(1, 0), state_plus(2, 0), state_plus(3, 0);
  //      quat_diff = quat_multiply(quat_st, Inv(quat_gt));
  //      double err_ori = (180 / M_PI) * 2 * quat_diff.block(0, 0, 3, 1).norm();
  //      // Calculate NEES values
  //      Eigen::Vector3d quat_diff_vec = quat_diff.block(0, 0, 3, 1);
  //      Eigen::Vector3d cov_vec = cov_plus.block(0, 0, 3, 3).inverse() * 2 * quat_diff.block(0, 0, 3, 1);
  //      double ori_nees = 2 * quat_diff_vec.dot(cov_vec);
  //      Eigen::Vector3d errpos = state_plus.block(4, 0, 3, 1) - state_gt.block(5, 0, 3, 1);
  //      double pos_nees = errpos.transpose() * cov_plus.block(3, 3, 3, 3).inverse() * errpos;
  //      PRINT_INFO(REDPURPLE "error to gt => %.3f, %.3f (deg,m) | nees => %.1f, %.1f (ori,pos) \n" RESET, err_ori, err_pos, ori_nees,
  //                 pos_nees);
  //    }
  //  }

  // Publish our odometry message if requested
  if (pub_odomimu.getNumSubscribers() != 0) {

    nav_msgs::Odometry odomIinM;
    odomIinM.header.stamp = ros::Time(timestamp);
    odomIinM.header.frame_id = "global";

    // The POSE component (orientation and position)
    odomIinM.pose.pose.orientation.x = state_plus(0);
    odomIinM.pose.pose.orientation.y = state_plus(1);
    odomIinM.pose.pose.orientation.z = state_plus(2);
    odomIinM.pose.pose.orientation.w = state_plus(3);
    odomIinM.pose.pose.position.x = state_plus(4);
    odomIinM.pose.pose.position.y = state_plus(5);
    odomIinM.pose.pose.position.z = state_plus(6);

    // The TWIST component (angular and linear velocities)
    odomIinM.child_frame_id = "imu";
    odomIinM.twist.twist.linear.x = state_plus(7);   // vel in local frame
    odomIinM.twist.twist.linear.y = state_plus(8);   // vel in local frame
    odomIinM.twist.twist.linear.z = state_plus(9);   // vel in local frame
    odomIinM.twist.twist.angular.x = state_plus(10); // we do not estimate this...
    odomIinM.twist.twist.angular.y = state_plus(11); // we do not estimate this...
    odomIinM.twist.twist.angular.z = state_plus(12); // we do not estimate this...

    // Finally set the covariance in the message (in the order position then orientation as per ros convention)
    Eigen::Matrix<double, 12, 12> Phi = Eigen::Matrix<double, 12, 12>::Zero();
    Phi.block(0, 3, 3, 3).setIdentity();
    Phi.block(3, 0, 3, 3).setIdentity();
    Phi.block(6, 6, 6, 6).setIdentity();
    cov_plus = Phi * cov_plus * Phi.transpose();
    for (int r = 0; r < 6; r++) {
      for (int c = 0; c < 6; c++) {
        odomIinM.pose.covariance[6 * r + c] = cov_plus(r, c);
      }
    }
    for (int r = 0; r < 6; r++) {
      for (int c = 0; c < 6; c++) {
        odomIinM.twist.covariance[6 * r + c] = cov_plus(r + 6, c + 6);
      }
    }
    pub_odomimu.publish(odomIinM);
  }

  // Publish our transform on TF
  // NOTE: since we use JPL we have an implicit conversion to Hamilton when we publish
  // NOTE: a rotation from GtoI in JPL has the same xyzw as a ItoG Hamilton rotation
  auto odom_pose = std::make_shared<ov_type::PoseJPL>();
  odom_pose->set_value(state_plus.block(0, 0, 7, 1));
  tf::StampedTransform trans = ROSVisualizerHelper::get_stamped_transform_from_pose(odom_pose, false);
  trans.frame_id_ = "global";
  trans.child_frame_id_ = "imu";
  if (publish_global2imu_tf) {
    mTfBr->sendTransform(trans);
  }

  // Loop through each camera calibration and publish it
  for (const auto &calib : state->_calib_IMUtoCAM) {
    tf::StampedTransform trans_calib = ROSVisualizerHelper::get_stamped_transform_from_pose(calib.second, true);
    trans_calib.frame_id_ = "imu";
    trans_calib.child_frame_id_ = "cam" + std::to_string(calib.first);
    if (publish_calibration_tf) {
      mTfBr->sendTransform(trans_calib);
    }
  }
}

// -----------------------------------------------------------------------------
// Final report: print calibration, RMSE/NEES, and elapsed time
// -----------------------------------------------------------------------------
void ROS1Visualizer::visualize_final() {

  // Final time offset value
  if (_app->get_state()->_options.do_calib_camera_timeoffset) {
    PRINT_INFO(REDPURPLE "camera-imu timeoffset = %.5f\n\n" RESET, _app->get_state()->_calib_dt_CAMtoIMU->value()(0));
  }

  // Final camera intrinsics
  if (_app->get_state()->_options.do_calib_camera_intrinsics) {
    for (int i = 0; i < _app->get_state()->_options.num_cameras; i++) {
      std::shared_ptr<Vec> calib = _app->get_state()->_cam_intrinsics.at(i);
      PRINT_INFO(REDPURPLE "cam%d intrinsics:\n" RESET, (int)i);
      PRINT_INFO(REDPURPLE "%.3f,%.3f,%.3f,%.3f\n" RESET, calib->value()(0), calib->value()(1), calib->value()(2), calib->value()(3));
      PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f,%.5f\n\n" RESET, calib->value()(4), calib->value()(5), calib->value()(6), calib->value()(7));
    }
  }

  // Final camera extrinsics
  if (_app->get_state()->_options.do_calib_camera_pose) {
    for (int i = 0; i < _app->get_state()->_options.num_cameras; i++) {
      std::shared_ptr<PoseJPL> calib = _app->get_state()->_calib_IMUtoCAM.at(i);
      Eigen::Matrix4d T_CtoI = Eigen::Matrix4d::Identity();
      T_CtoI.block(0, 0, 3, 3) = quat_2_Rot(calib->quat()).transpose();
      T_CtoI.block(0, 3, 3, 1) = -T_CtoI.block(0, 0, 3, 3) * calib->pos();
      PRINT_INFO(REDPURPLE "T_C%dtoI:\n" RESET, i);
      PRINT_INFO(REDPURPLE "%.3f,%.3f,%.3f,%.3f,\n" RESET, T_CtoI(0, 0), T_CtoI(0, 1), T_CtoI(0, 2), T_CtoI(0, 3));
      PRINT_INFO(REDPURPLE "%.3f,%.3f,%.3f,%.3f,\n" RESET, T_CtoI(1, 0), T_CtoI(1, 1), T_CtoI(1, 2), T_CtoI(1, 3));
      PRINT_INFO(REDPURPLE "%.3f,%.3f,%.3f,%.3f,\n" RESET, T_CtoI(2, 0), T_CtoI(2, 1), T_CtoI(2, 2), T_CtoI(2, 3));
      PRINT_INFO(REDPURPLE "%.3f,%.3f,%.3f,%.3f\n\n" RESET, T_CtoI(3, 0), T_CtoI(3, 1), T_CtoI(3, 2), T_CtoI(3, 3));
    }
  }

  // IMU intrinsics
  if (_app->get_state()->_options.do_calib_imu_intrinsics) {
    Eigen::Matrix3d Dw = State::Dm(_app->get_state()->_options.imu_model, _app->get_state()->_calib_imu_dw->value());
    Eigen::Matrix3d Da = State::Dm(_app->get_state()->_options.imu_model, _app->get_state()->_calib_imu_da->value());
    Eigen::Matrix3d Tw = Dw.colPivHouseholderQr().solve(Eigen::Matrix3d::Identity());
    Eigen::Matrix3d Ta = Da.colPivHouseholderQr().solve(Eigen::Matrix3d::Identity());
    Eigen::Matrix3d R_IMUtoACC = _app->get_state()->_calib_imu_ACCtoIMU->Rot().transpose();
    Eigen::Matrix3d R_IMUtoGYRO = _app->get_state()->_calib_imu_GYROtoIMU->Rot().transpose();
    PRINT_INFO(REDPURPLE "Tw:\n" RESET);
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f,\n" RESET, Tw(0, 0), Tw(0, 1), Tw(0, 2));
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f,\n" RESET, Tw(1, 0), Tw(1, 1), Tw(1, 2));
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f\n\n" RESET, Tw(2, 0), Tw(2, 1), Tw(2, 2));
    PRINT_INFO(REDPURPLE "Ta:\n" RESET);
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f,\n" RESET, Ta(0, 0), Ta(0, 1), Ta(0, 2));
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f,\n" RESET, Ta(1, 0), Ta(1, 1), Ta(1, 2));
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f\n\n" RESET, Ta(2, 0), Ta(2, 1), Ta(2, 2));
    PRINT_INFO(REDPURPLE "R_IMUtoACC:\n" RESET);
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f,\n" RESET, R_IMUtoACC(0, 0), R_IMUtoACC(0, 1), R_IMUtoACC(0, 2));
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f,\n" RESET, R_IMUtoACC(1, 0), R_IMUtoACC(1, 1), R_IMUtoACC(1, 2));
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f\n\n" RESET, R_IMUtoACC(2, 0), R_IMUtoACC(2, 1), R_IMUtoACC(2, 2));
    PRINT_INFO(REDPURPLE "R_IMUtoGYRO:\n" RESET);
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f,\n" RESET, R_IMUtoGYRO(0, 0), R_IMUtoGYRO(0, 1), R_IMUtoGYRO(0, 2));
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f,\n" RESET, R_IMUtoGYRO(1, 0), R_IMUtoGYRO(1, 1), R_IMUtoGYRO(1, 2));
    PRINT_INFO(REDPURPLE "%.5f,%.5f,%.5f\n\n" RESET, R_IMUtoGYRO(2, 0), R_IMUtoGYRO(2, 1), R_IMUtoGYRO(2, 2));
  }

  // IMU intrinsics gravity sensitivity
  if (_app->get_state()->_options.do_calib_imu_g_sensitivity) {
    Eigen::Matrix3d Tg = State::Tg(_app->get_state()->_calib_imu_tg->value());
    PRINT_INFO(REDPURPLE "Tg:\n" RESET);
    PRINT_INFO(REDPURPLE "%.6f,%.6f,%.6f,\n" RESET, Tg(0, 0), Tg(0, 1), Tg(0, 2));
    PRINT_INFO(REDPURPLE "%.6f,%.6f,%.6f,\n" RESET, Tg(1, 0), Tg(1, 1), Tg(1, 2));
    PRINT_INFO(REDPURPLE "%.6f,%.6f,%.6f\n\n" RESET, Tg(2, 0), Tg(2, 1), Tg(2, 2));
  }

  // Publish RMSE if we have it
  if (!gt_states.empty()) {
    PRINT_INFO(REDPURPLE "RMSE: %.3f (deg) orientation\n" RESET, std::sqrt(summed_mse_ori / summed_number));
    PRINT_INFO(REDPURPLE "RMSE: %.3f (m) position\n\n" RESET, std::sqrt(summed_mse_pos / summed_number));
  }

  // Publish RMSE and NEES if doing simulation
  if (_sim != nullptr) {
    PRINT_INFO(REDPURPLE "RMSE: %.3f (deg) orientation\n" RESET, std::sqrt(summed_mse_ori / summed_number));
    PRINT_INFO(REDPURPLE "RMSE: %.3f (m) position\n\n" RESET, std::sqrt(summed_mse_pos / summed_number));
    PRINT_INFO(REDPURPLE "NEES: %.3f (deg) orientation\n" RESET, summed_nees_ori / summed_number);
    PRINT_INFO(REDPURPLE "NEES: %.3f (m) position\n\n" RESET, summed_nees_pos / summed_number);
  }

  // Print the total time
  rT2 = boost::posix_time::microsec_clock::local_time();
  PRINT_INFO(REDPURPLE "TIME: %.3f seconds\n\n" RESET, (rT2 - rT1).total_microseconds() * 1e-6);
}

// -----------------------------------------------------------------------------
// Input callback: IMU propagation and queued camera update release
// -----------------------------------------------------------------------------
void ROS1Visualizer::callback_inertial(const sensor_msgs::Imu::ConstPtr &msg) {

  // ROS sensor_msgs/Imu 메시지를 OpenVINS 내부 형식인 ov_core::ImuData로 변환한다.
  // timestamp는 ROS header stamp를 초 단위 double로 바꾼 값이다.
  ov_core::ImuData message;
  message.timestamp = msg->header.stamp.toSec();

  // wm: angular velocity [rad/s], am: linear acceleration [m/s^2].
  // ROS 메시지의 x/y/z 값을 Eigen vector 형태로 복사한다.
  message.wm << msg->angular_velocity.x, msg->angular_velocity.y, msg->angular_velocity.z;
  message.am << msg->linear_acceleration.x, msg->linear_acceleration.y, msg->linear_acceleration.z;

  // IMU 측정값을 VIO estimator에 넣어 propagation을 수행한다.
  // 카메라 update 사이에서도 IMU가 들어올 때마다 상태가 계속 전파된다.
  _app->feed_measurement_imu(message);

  // 최신 IMU timestamp 기준으로 high-rate odometry를 publish/visualize한다.
  // 카메라 update보다 높은 주기로 IMU propagation 결과를 볼 수 있게 해준다.
  visualize_odometry(message.timestamp);

  // 이미 다른 thread가 camera_queue를 처리 중이면, 이번 IMU callback에서는 추가 작업을 하지 않는다.
  // 이렇게 해야 IMU callback이 오래 막히지 않고 다음 IMU 측정을 계속 받을 수 있다.
  if (thread_update_running)
    return;

  // 여기부터는 아직 camera update 처리 thread가 없다는 뜻이다.
  // flag를 세워 중복 thread 생성을 막고, 아래에서 큐 처리 thread를 만든다.
  thread_update_running = true;
  std::thread thread([&] {
    // camera_queue를 처리하는 동안 새 image callback이 queue를 동시에 수정하지 못하도록 잠근다.
    // 이 lock 범위 안에서 queue 검사, camera update, pop_front가 모두 일어난다.
    std::lock_guard<std::mutex> lck(camera_queue_mtx);

    // 현재 queue 안에 들어와 있는 image stream 종류를 센다.
    // sensor_ids.at(0)는 해당 CameraData가 어느 카메라 stream에서 왔는지 나타낸다.
    std::map<int, bool> unique_cam_ids;
    for (const auto &cam_msg : camera_queue) {
      unique_cam_ids[cam_msg.sensor_ids.at(0)] = true;
    }

    // 필요한 카메라 stream이 queue에 모두 들어올 때까지 기다린다.
    // stereo 설정(num_cameras == 2)에서는 좌/우 이미지가 하나의 CameraData로 묶여 들어오므로
    // unique stream 개수를 1개로 본다. 그 외에는 설정된 카메라 수만큼 기다린다.
    auto params = _app->get_params();
    size_t num_unique_cameras = (params.state_options.num_cameras == 2) ? 1 : params.state_options.num_cameras;
    if (unique_cam_ids.size() == num_unique_cameras) {

      // 현재 IMU timestamp를 camera clock 기준으로 변환한다.
      // OpenVINS는 camera timestamp의 측정을 처리하려면 그 시각 이후의 IMU 측정이 적어도 하나 필요하다.
      double timestamp_imu_inC = message.timestamp - _app->get_state()->_calib_dt_CAMtoIMU->value()(0);

      // queue의 가장 오래된 camera measurement부터 처리 가능한지 확인한다.
      // camera timestamp가 현재 IMU-in-camera-clock보다 과거이면, 충분한 IMU propagation이 가능하므로 update한다.
      while (!camera_queue.empty() && camera_queue.at(0).timestamp < timestamp_imu_inC) {
        auto rT0_1 = boost::posix_time::microsec_clock::local_time();

        // update_dt는 처리 지연 정도를 ms 단위로 출력하기 위한 값이다.
        // 아래 PRINT_INFO에서는 "%.2f ms behind"로 표시된다.
        double update_dt = 100.0 * (timestamp_imu_inC - camera_queue.at(0).timestamp);

        // camera measurement를 estimator에 넣어 MSCKF/SLAM update를 수행한다.
        // 이 호출 이후 state는 카메라 정보가 반영된 posterior 상태가 된다.
        _app->feed_measurement_camera(camera_queue.at(0));

        // update된 state, feature, image 등을 ROS topic으로 publish한다.
        visualize();

        // 처리가 끝난 camera measurement는 queue에서 제거한다.
        camera_queue.pop_front();

        // 이번 camera update와 visualization에 걸린 시간을 계산해 로그로 출력한다.
        auto rT0_2 = boost::posix_time::microsec_clock::local_time();
        double time_total = (rT0_2 - rT0_1).total_microseconds() * 1e-6;
        PRINT_INFO(BLUE "[TIME]: %.4f seconds total (%.1f hz, %.2f ms behind)\n" RESET, time_total, 1.0 / time_total, update_dt);
      }
    }

    // queue 처리 thread가 끝났으므로 다음 IMU callback에서 새 update thread를 만들 수 있게 한다.
    thread_update_running = false;
  });

  // 설정이 single-threaded이면 여기서 thread가 끝날 때까지 기다린다.
  // multi-threaded이면 detach해서 background에서 camera update가 돌게 하고,
  // 현재 IMU callback은 바로 반환되어 다음 ROS callback을 받을 수 있게 한다.
  if (!_app->get_params().use_multi_threading_subs) {
    thread.join();
  } else {
    thread.detach();
  }
}

// -----------------------------------------------------------------------------
// Input callback: monocular image conversion and camera queue insertion
// -----------------------------------------------------------------------------
void ROS1Visualizer::callback_monocular(const sensor_msgs::ImageConstPtr &msg0, int cam_id0) {

  // Check if we should drop this image
  double timestamp = msg0->header.stamp.toSec();
  double time_delta = 1.0 / _app->get_params().track_frequency;
  if (camera_last_timestamp.find(cam_id0) != camera_last_timestamp.end() && timestamp < camera_last_timestamp.at(cam_id0) + time_delta) {
    return;
  }
  camera_last_timestamp[cam_id0] = timestamp;

  // Get the image
  cv_bridge::CvImageConstPtr cv_ptr;
  try {
    cv_ptr = cv_bridge::toCvShare(msg0, sensor_msgs::image_encodings::MONO8);
  } catch (cv_bridge::Exception &e) {
    PRINT_ERROR("cv_bridge exception: %s", e.what());
    return;
  }

  // Create the measurement
  ov_core::CameraData message;
  message.timestamp = cv_ptr->header.stamp.toSec();
  message.sensor_ids.push_back(cam_id0);
  message.images.push_back(cv_ptr->image.clone());

  // Load the mask if we are using it, else it is empty
  // TODO: in the future we should get this from external pixel segmentation
  if (_app->get_params().use_mask) {
    message.masks.push_back(_app->get_params().masks.at(cam_id0));
  } else {
    message.masks.push_back(cv::Mat::zeros(cv_ptr->image.rows, cv_ptr->image.cols, CV_8UC1));
  }

  // append it to our queue of images
  std::lock_guard<std::mutex> lck(camera_queue_mtx);
  camera_queue.push_back(message);
  std::sort(camera_queue.begin(), camera_queue.end());
}

// -----------------------------------------------------------------------------
// Input callback: stereo image conversion and camera queue insertion
// -----------------------------------------------------------------------------
void ROS1Visualizer::callback_stereo(const sensor_msgs::ImageConstPtr &msg0, const sensor_msgs::ImageConstPtr &msg1, int cam_id0,
                                     int cam_id1) {

  // Check if we should drop this image
  double timestamp = msg0->header.stamp.toSec();
  double time_delta = 1.0 / _app->get_params().track_frequency;
  if (camera_last_timestamp.find(cam_id0) != camera_last_timestamp.end() && timestamp < camera_last_timestamp.at(cam_id0) + time_delta) {
    return;
  }
  camera_last_timestamp[cam_id0] = timestamp;

  // Get the image
  cv_bridge::CvImageConstPtr cv_ptr0;
  try {
    cv_ptr0 = cv_bridge::toCvShare(msg0, sensor_msgs::image_encodings::MONO8);
  } catch (cv_bridge::Exception &e) {
    PRINT_ERROR("cv_bridge exception: %s\n", e.what());
    return;
  }

  // Get the image
  cv_bridge::CvImageConstPtr cv_ptr1;
  try {
    cv_ptr1 = cv_bridge::toCvShare(msg1, sensor_msgs::image_encodings::MONO8);
  } catch (cv_bridge::Exception &e) {
    PRINT_ERROR("cv_bridge exception: %s\n", e.what());
    return;
  }

  // Create the measurement
  ov_core::CameraData message;
  message.timestamp = cv_ptr0->header.stamp.toSec();
  message.sensor_ids.push_back(cam_id0);
  message.sensor_ids.push_back(cam_id1);
  message.images.push_back(cv_ptr0->image.clone());
  message.images.push_back(cv_ptr1->image.clone());

  // Load the mask if we are using it, else it is empty
  // TODO: in the future we should get this from external pixel segmentation
  if (_app->get_params().use_mask) {
    message.masks.push_back(_app->get_params().masks.at(cam_id0));
    message.masks.push_back(_app->get_params().masks.at(cam_id1));
  } else {
    // message.masks.push_back(cv::Mat(cv_ptr0->image.rows, cv_ptr0->image.cols, CV_8UC1, cv::Scalar(255)));
    message.masks.push_back(cv::Mat::zeros(cv_ptr0->image.rows, cv_ptr0->image.cols, CV_8UC1));
    message.masks.push_back(cv::Mat::zeros(cv_ptr1->image.rows, cv_ptr1->image.cols, CV_8UC1));
  }

  // append it to our queue of images
  std::lock_guard<std::mutex> lck(camera_queue_mtx);
  camera_queue.push_back(message);
  std::sort(camera_queue.begin(), camera_queue.end());
}

// -----------------------------------------------------------------------------
// Output publisher: current IMU pose with covariance and accumulated path
// -----------------------------------------------------------------------------
void ROS1Visualizer::publish_state() {

  // Get the current state
  std::shared_ptr<State> state = _app->get_state();

  // We want to publish in the IMU clock frame
  // The timestamp in the state will be the last camera time
  double t_ItoC = state->_calib_dt_CAMtoIMU->value()(0);
  double timestamp_inI = state->_timestamp + t_ItoC;

  // Create pose of IMU (note we use the bag time)
  geometry_msgs::PoseWithCovarianceStamped poseIinM;
  poseIinM.header.stamp = ros::Time(timestamp_inI);
  poseIinM.header.seq = poses_seq_imu;
  poseIinM.header.frame_id = "global";
  poseIinM.pose.pose.orientation.x = state->_imu->quat()(0);
  poseIinM.pose.pose.orientation.y = state->_imu->quat()(1);
  poseIinM.pose.pose.orientation.z = state->_imu->quat()(2);
  poseIinM.pose.pose.orientation.w = state->_imu->quat()(3);
  poseIinM.pose.pose.position.x = state->_imu->pos()(0);
  poseIinM.pose.pose.position.y = state->_imu->pos()(1);
  poseIinM.pose.pose.position.z = state->_imu->pos()(2);

  // Finally set the covariance in the message (in the order position then orientation as per ros convention)
  std::vector<std::shared_ptr<Type>> statevars;
  statevars.push_back(state->_imu->pose()->p());
  statevars.push_back(state->_imu->pose()->q());
  Eigen::Matrix<double, 6, 6> covariance_posori = StateHelper::get_marginal_covariance(_app->get_state(), statevars);
  for (int r = 0; r < 6; r++) {
    for (int c = 0; c < 6; c++) {
      poseIinM.pose.covariance[6 * r + c] = covariance_posori(r, c);
    }
  }
  pub_poseimu.publish(poseIinM);

  //=========================================================
  //=========================================================

  // Append to our pose vector
  geometry_msgs::PoseStamped posetemp;
  posetemp.header = poseIinM.header;
  posetemp.pose = poseIinM.pose.pose;
  poses_imu.push_back(posetemp);

  // Create our path (imu)
  // NOTE: We downsample the number of poses as needed to prevent rviz crashes
  // NOTE: https://github.com/ros-visualization/rviz/issues/1107
  nav_msgs::Path arrIMU;
  arrIMU.header.stamp = ros::Time::now();
  arrIMU.header.seq = poses_seq_imu;
  arrIMU.header.frame_id = "global";
  for (size_t i = 0; i < poses_imu.size(); i += std::floor((double)poses_imu.size() / 16384.0) + 1) {
    arrIMU.poses.push_back(poses_imu.at(i));
  }
  pub_pathimu.publish(arrIMU);

  // Move them forward in time
  poses_seq_imu++;
}

// -----------------------------------------------------------------------------
// Output publisher: tracker history image
// -----------------------------------------------------------------------------
void ROS1Visualizer::publish_images() {

  // Return if we have already visualized
  if (_app->get_state() == nullptr)
    return;
  if (last_visualization_timestamp_image == _app->get_state()->_timestamp && _app->initialized())
    return;
  last_visualization_timestamp_image = _app->get_state()->_timestamp;

  // Check if we have subscribers
  if (it_pub_tracks.getNumSubscribers() == 0)
    return;

  // Get our image of history tracks
  cv::Mat img_history = _app->get_historical_viz_image();
  if (img_history.empty())
    return;

  // Create our message
  std_msgs::Header header;
  header.stamp = ros::Time::now();
  header.frame_id = "cam0";
  sensor_msgs::ImagePtr exl_msg = cv_bridge::CvImage(header, "bgr8", img_history).toImageMsg();

  // Publish
  it_pub_tracks.publish(exl_msg);
}

// -----------------------------------------------------------------------------
// Output publisher: MSCKF, SLAM, ARUCO, and simulation feature clouds
// -----------------------------------------------------------------------------
void ROS1Visualizer::publish_features() {

  // Check if we have subscribers
  if (pub_points_msckf.getNumSubscribers() == 0 && pub_points_slam.getNumSubscribers() == 0 && pub_points_aruco.getNumSubscribers() == 0 &&
      pub_points_sim.getNumSubscribers() == 0)
    return;

  // Get our good MSCKF features
  std::vector<Eigen::Vector3d> feats_msckf = _app->get_good_features_MSCKF();
  sensor_msgs::PointCloud2 cloud = ROSVisualizerHelper::get_ros_pointcloud(feats_msckf);
  pub_points_msckf.publish(cloud);

  // Get our good SLAM features
  std::vector<Eigen::Vector3d> feats_slam = _app->get_features_SLAM();
  sensor_msgs::PointCloud2 cloud_SLAM = ROSVisualizerHelper::get_ros_pointcloud(feats_slam);
  pub_points_slam.publish(cloud_SLAM);

  // Get our good ARUCO features
  std::vector<Eigen::Vector3d> feats_aruco = _app->get_features_ARUCO();
  sensor_msgs::PointCloud2 cloud_ARUCO = ROSVisualizerHelper::get_ros_pointcloud(feats_aruco);
  pub_points_aruco.publish(cloud_ARUCO);

  // Skip the rest of we are not doing simulation
  if (_sim == nullptr)
    return;

  // Get our good SIMULATION features
  std::vector<Eigen::Vector3d> feats_sim = _sim->get_map_vec();
  sensor_msgs::PointCloud2 cloud_SIM = ROSVisualizerHelper::get_ros_pointcloud(feats_sim);
  pub_points_sim.publish(cloud_SIM);
}

// -----------------------------------------------------------------------------
// Output publisher: groundtruth pose/path plus RMSE and NEES bookkeeping
// -----------------------------------------------------------------------------
void ROS1Visualizer::publish_groundtruth() {

  // Our groundtruth state
  Eigen::Matrix<double, 17, 1> state_gt;

  // We want to publish in the IMU clock frame
  // The timestamp in the state will be the last camera time
  double t_ItoC = _app->get_state()->_calib_dt_CAMtoIMU->value()(0);
  double timestamp_inI = _app->get_state()->_timestamp + t_ItoC;

  // Check that we have the timestamp in our GT file [time(sec),q_GtoI,p_IinG,v_IinG,b_gyro,b_accel]
  if (_sim == nullptr && (gt_states.empty() || !DatasetReader::get_gt_state(timestamp_inI, state_gt, gt_states))) {
    return;
  }

  // Get the simulated groundtruth
  // NOTE: we get the true time in the IMU clock frame
  if (_sim != nullptr) {
    timestamp_inI = _app->get_state()->_timestamp + _sim->get_true_parameters().calib_camimu_dt;
    if (!_sim->get_state(timestamp_inI, state_gt))
      return;
  }

  // Get the GT and system state state
  Eigen::Matrix<double, 16, 1> state_ekf = _app->get_state()->_imu->value();

  // Create pose of IMU
  geometry_msgs::PoseStamped poseIinM;
  poseIinM.header.stamp = ros::Time(timestamp_inI);
  poseIinM.header.seq = poses_seq_gt;
  poseIinM.header.frame_id = "global";
  poseIinM.pose.orientation.x = state_gt(1, 0);
  poseIinM.pose.orientation.y = state_gt(2, 0);
  poseIinM.pose.orientation.z = state_gt(3, 0);
  poseIinM.pose.orientation.w = state_gt(4, 0);
  poseIinM.pose.position.x = state_gt(5, 0);
  poseIinM.pose.position.y = state_gt(6, 0);
  poseIinM.pose.position.z = state_gt(7, 0);
  pub_posegt.publish(poseIinM);

  // Append to our pose vector
  poses_gt.push_back(poseIinM);

  // Create our path (imu)
  // NOTE: We downsample the number of poses as needed to prevent rviz crashes
  // NOTE: https://github.com/ros-visualization/rviz/issues/1107
  nav_msgs::Path arrIMU;
  arrIMU.header.stamp = ros::Time::now();
  arrIMU.header.seq = poses_seq_gt;
  arrIMU.header.frame_id = "global";
  for (size_t i = 0; i < poses_gt.size(); i += std::floor((double)poses_gt.size() / 16384.0) + 1) {
    arrIMU.poses.push_back(poses_gt.at(i));
  }
  pub_pathgt.publish(arrIMU);

  // Move them forward in time
  poses_seq_gt++;

  // Publish our transform on TF
  tf::StampedTransform trans;
  trans.stamp_ = ros::Time::now();
  trans.frame_id_ = "global";
  trans.child_frame_id_ = "truth";
  tf::Quaternion quat(state_gt(1, 0), state_gt(2, 0), state_gt(3, 0), state_gt(4, 0));
  trans.setRotation(quat);
  tf::Vector3 orig(state_gt(5, 0), state_gt(6, 0), state_gt(7, 0));
  trans.setOrigin(orig);
  if (publish_global2imu_tf) {
    mTfBr->sendTransform(trans);
  }

  //==========================================================================
  //==========================================================================

  // Difference between positions
  double dx = state_ekf(4, 0) - state_gt(5, 0);
  double dy = state_ekf(5, 0) - state_gt(6, 0);
  double dz = state_ekf(6, 0) - state_gt(7, 0);
  double err_pos = std::sqrt(dx * dx + dy * dy + dz * dz);

  // Quaternion error
  Eigen::Matrix<double, 4, 1> quat_gt, quat_st, quat_diff;
  quat_gt << state_gt(1, 0), state_gt(2, 0), state_gt(3, 0), state_gt(4, 0);
  quat_st << state_ekf(0, 0), state_ekf(1, 0), state_ekf(2, 0), state_ekf(3, 0);
  quat_diff = quat_multiply(quat_st, Inv(quat_gt));
  double err_ori = (180 / M_PI) * 2 * quat_diff.block(0, 0, 3, 1).norm();

  //==========================================================================
  //==========================================================================

  // Get covariance of pose
  std::vector<std::shared_ptr<Type>> statevars;
  statevars.push_back(_app->get_state()->_imu->q());
  statevars.push_back(_app->get_state()->_imu->p());
  Eigen::Matrix<double, 6, 6> covariance = StateHelper::get_marginal_covariance(_app->get_state(), statevars);

  // Calculate NEES values
  // NOTE: need to manually multiply things out to make static asserts work
  // NOTE: https://github.com/rpng/open_vins/pull/226
  // NOTE: https://github.com/rpng/open_vins/issues/236
  // NOTE: https://gitlab.com/libeigen/eigen/-/issues/1664
  Eigen::Vector3d quat_diff_vec = quat_diff.block(0, 0, 3, 1);
  Eigen::Vector3d cov_vec = covariance.block(0, 0, 3, 3).inverse() * 2 * quat_diff.block(0, 0, 3, 1);
  double ori_nees = 2 * quat_diff_vec.dot(cov_vec);
  Eigen::Vector3d errpos = state_ekf.block(4, 0, 3, 1) - state_gt.block(5, 0, 3, 1);
  double pos_nees = errpos.transpose() * covariance.block(3, 3, 3, 3).inverse() * errpos;

  //==========================================================================
  //==========================================================================

  // Update our average variables
  if (!std::isnan(ori_nees) && !std::isnan(pos_nees)) {
    summed_mse_ori += err_ori * err_ori;
    summed_mse_pos += err_pos * err_pos;
    summed_nees_ori += ori_nees;
    summed_nees_pos += pos_nees;
    summed_number++;
  }

  // Nice display for the user
  PRINT_INFO(REDPURPLE "error to gt => %.3f, %.3f (deg,m) | rmse => %.3f, %.3f (deg,m) | called %d times\n" RESET, err_ori, err_pos,
             std::sqrt(summed_mse_ori / summed_number), std::sqrt(summed_mse_pos / summed_number), (int)summed_number);
  PRINT_INFO(REDPURPLE "nees => %.1f, %.1f (ori,pos) | avg nees = %.1f, %.1f (ori,pos)\n" RESET, ori_nees, pos_nees,
             summed_nees_ori / summed_number, summed_nees_pos / summed_number);

  //==========================================================================
  //==========================================================================
}

// -----------------------------------------------------------------------------
// Output publisher: loop-closure pose, calibration, feature, and depth products
// -----------------------------------------------------------------------------
void ROS1Visualizer::publish_loopclosure_information() {

  // Get the current tracks in this frame
  double active_tracks_time1 = -1;
  double active_tracks_time2 = -1;
  std::unordered_map<size_t, Eigen::Vector3d> active_tracks_posinG;
  std::unordered_map<size_t, Eigen::Vector3d> active_tracks_uvd;
  cv::Mat active_cam0_image;
  _app->get_active_tracks(active_tracks_time1, active_tracks_posinG, active_tracks_uvd);
  _app->get_active_image(active_tracks_time2, active_cam0_image);
  if (active_tracks_time1 == -1)
    return;
  if (_app->get_state()->_clones_IMU.find(active_tracks_time1) == _app->get_state()->_clones_IMU.end())
    return;
  Eigen::Vector4d quat = _app->get_state()->_clones_IMU.at(active_tracks_time1)->quat();
  Eigen::Vector3d pos = _app->get_state()->_clones_IMU.at(active_tracks_time1)->pos();
  if (active_tracks_time1 != active_tracks_time2)
    return;

  // Default header
  std_msgs::Header header;
  header.stamp = ros::Time(active_tracks_time1);

  //======================================================
  // Check if we have subscribers for the pose odometry, camera intrinsics, or extrinsics
  if (pub_loop_pose.getNumSubscribers() != 0 || pub_loop_extrinsic.getNumSubscribers() != 0 ||
      pub_loop_intrinsics.getNumSubscribers() != 0) {

    // PUBLISH HISTORICAL POSE ESTIMATE
    nav_msgs::Odometry odometry_pose;
    odometry_pose.header = header;
    odometry_pose.header.frame_id = "global";
    odometry_pose.pose.pose.position.x = pos(0);
    odometry_pose.pose.pose.position.y = pos(1);
    odometry_pose.pose.pose.position.z = pos(2);
    odometry_pose.pose.pose.orientation.x = quat(0);
    odometry_pose.pose.pose.orientation.y = quat(1);
    odometry_pose.pose.pose.orientation.z = quat(2);
    odometry_pose.pose.pose.orientation.w = quat(3);
    pub_loop_pose.publish(odometry_pose);

    // PUBLISH IMU TO CAMERA0 EXTRINSIC
    // need to flip the transform to the IMU frame
    Eigen::Vector4d q_ItoC = _app->get_state()->_calib_IMUtoCAM.at(0)->quat();
    Eigen::Vector3d p_CinI = -_app->get_state()->_calib_IMUtoCAM.at(0)->Rot().transpose() * _app->get_state()->_calib_IMUtoCAM.at(0)->pos();
    nav_msgs::Odometry odometry_calib;
    odometry_calib.header = header;
    odometry_calib.header.frame_id = "imu";
    odometry_calib.pose.pose.position.x = p_CinI(0);
    odometry_calib.pose.pose.position.y = p_CinI(1);
    odometry_calib.pose.pose.position.z = p_CinI(2);
    odometry_calib.pose.pose.orientation.x = q_ItoC(0);
    odometry_calib.pose.pose.orientation.y = q_ItoC(1);
    odometry_calib.pose.pose.orientation.z = q_ItoC(2);
    odometry_calib.pose.pose.orientation.w = q_ItoC(3);
    pub_loop_extrinsic.publish(odometry_calib);

    // PUBLISH CAMERA0 INTRINSICS
    bool is_fisheye = (std::dynamic_pointer_cast<ov_core::CamEqui>(_app->get_params().camera_intrinsics.at(0)) != nullptr);
    sensor_msgs::CameraInfo cameraparams;
    cameraparams.header = header;
    cameraparams.header.frame_id = "cam0";
    cameraparams.distortion_model = is_fisheye ? "equidistant" : "plumb_bob";
    Eigen::VectorXd cparams = _app->get_state()->_cam_intrinsics.at(0)->value();
    cameraparams.D = {cparams(4), cparams(5), cparams(6), cparams(7)};
    cameraparams.K = {cparams(0), 0, cparams(2), 0, cparams(1), cparams(3), 0, 0, 1};
    pub_loop_intrinsics.publish(cameraparams);
  }

  //======================================================
  // PUBLISH FEATURE TRACKS IN THE GLOBAL FRAME OF REFERENCE
  if (pub_loop_point.getNumSubscribers() != 0) {

    // Construct the message
    sensor_msgs::PointCloud point_cloud;
    point_cloud.header = header;
    point_cloud.header.frame_id = "global";
    for (const auto &feattimes : active_tracks_posinG) {

      // Get this feature information
      size_t featid = feattimes.first;
      Eigen::Vector3d uvd = Eigen::Vector3d::Zero();
      if (active_tracks_uvd.find(featid) != active_tracks_uvd.end()) {
        uvd = active_tracks_uvd.at(featid);
      }
      Eigen::Vector3d pFinG = active_tracks_posinG.at(featid);

      // Push back 3d point
      geometry_msgs::Point32 p;
      p.x = pFinG(0);
      p.y = pFinG(1);
      p.z = pFinG(2);
      point_cloud.points.push_back(p);

      // Push back the uv_norm, uv_raw, and feature id
      // NOTE: we don't use the normalized coordinates to save time here
      // NOTE: they will have to be re-normalized in the loop closure code
      sensor_msgs::ChannelFloat32 p_2d;
      p_2d.values.push_back(0);
      p_2d.values.push_back(0);
      p_2d.values.push_back(uvd(0));
      p_2d.values.push_back(uvd(1));
      p_2d.values.push_back(featid);
      point_cloud.channels.push_back(p_2d);
    }
    pub_loop_point.publish(point_cloud);
  }

  //======================================================
  // Depth images of sparse points and its colorized version
  if (it_pub_loop_img_depth.getNumSubscribers() != 0 || it_pub_loop_img_depth_color.getNumSubscribers() != 0) {

    // Create the images we will populate with the depths
    std::pair<int, int> wh_pair = {active_cam0_image.cols, active_cam0_image.rows};
    cv::Mat depthmap = cv::Mat::zeros(wh_pair.second, wh_pair.first, CV_16UC1);
    cv::Mat depthmap_viz = active_cam0_image;

    // Loop through all points and append
    for (const auto &feattimes : active_tracks_uvd) {

      // Get this feature information
      size_t featid = feattimes.first;
      Eigen::Vector3d uvd = active_tracks_uvd.at(featid);

      // Skip invalid points
      double dw = 4;
      if (uvd(0) < dw || uvd(0) > wh_pair.first - dw || uvd(1) < dw || uvd(1) > wh_pair.second - dw) {
        continue;
      }

      // Append the depth
      // NOTE: scaled by 1000 to fit the 16U
      // NOTE: access order is y,x (stupid opencv convention stuff)
      depthmap.at<uint16_t>((int)uvd(1), (int)uvd(0)) = (uint16_t)(1000 * uvd(2));

      // Taken from LSD-SLAM codebase segment into 0-4 meter segments:
      // https://github.com/tum-vision/lsd_slam/blob/d1e6f0e1a027889985d2e6b4c0fe7a90b0c75067/lsd_slam_core/src/util/globalFuncs.cpp#L87-L96
      float id = 1.0f / (float)uvd(2);
      float r = (0.0f - id) * 255 / 1.0f;
      if (r < 0)
        r = -r;
      float g = (1.0f - id) * 255 / 1.0f;
      if (g < 0)
        g = -g;
      float b = (2.0f - id) * 255 / 1.0f;
      if (b < 0)
        b = -b;
      uchar rc = r < 0 ? 0 : (r > 255 ? 255 : r);
      uchar gc = g < 0 ? 0 : (g > 255 ? 255 : g);
      uchar bc = b < 0 ? 0 : (b > 255 ? 255 : b);
      cv::Scalar color(255 - rc, 255 - gc, 255 - bc);

      // Small square around the point (note the above bound check needs to take into account this width)
      cv::Point p0(uvd(0) - dw, uvd(1) - dw);
      cv::Point p1(uvd(0) + dw, uvd(1) + dw);
      cv::rectangle(depthmap_viz, p0, p1, color, -1);
    }

    // Create our messages
    header.frame_id = "cam0";
    sensor_msgs::ImagePtr exl_msg1 = cv_bridge::CvImage(header, sensor_msgs::image_encodings::TYPE_16UC1, depthmap).toImageMsg();
    it_pub_loop_img_depth.publish(exl_msg1);
    header.stamp = ros::Time::now();
    header.frame_id = "cam0";
    sensor_msgs::ImagePtr exl_msg2 = cv_bridge::CvImage(header, "bgr8", depthmap_viz).toImageMsg();
    it_pub_loop_img_depth_color.publish(exl_msg2);
  }
}
