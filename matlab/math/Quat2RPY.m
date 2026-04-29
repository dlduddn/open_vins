function euler_ang = Quat2RPY(orient)
    % Quat2RPY: 쿼터니언을 roll/pitch/yaw 오일러 각으로 변환
    % 입력: orient - SO3 형태의 쿼터니언 (4x1 벡터: [qw, qx, qy, qz]')
    % 출력: euler_ang - 오일러 각 radian (3x1 벡터: [roll, pitch, yaw]')

    % 쿼터니언 값 추출
    qw = orient(1);
    qx = orient(2);
    qy = orient(3);
    qz = orient(4);

    % 사전 계산된 쿼터니언 값의 제곱
    sqw = qw * qw;
    sqx = qx * qx;
    sqy = qy * qy;
    sqz = qz * qz;

    % 쿼터니언 정규화 확인
    unit = sqx + sqy + sqz + sqw; % 정규화 단위
    test = qw * qy - qz * qx;    % 극점 여부 확인

    % 극점(특이점) 처리
    if test > 0.49999 * unit
        % 북극 특이점
        roll = 2 * atan2(qx, qw);
        pitch = pi / 2;
        yaw = 0;
        euler_ang = [roll, pitch, yaw]';
        return;
    elseif test < -0.49999 * unit
        % 남극 특이점
        roll = -2 * atan2(qx, qw);
        pitch = -pi / 2;
        yaw = 0;
        euler_ang = [roll, pitch, yaw]';
        return;
    end

    % 일반적인 경우의 오일러 각 계산
    roll = atan2(2 * (qx * qw + qy * qz), -sqx - sqy + sqz + sqw);
    pitch = asin(2 * test / unit);
    yaw = atan2(2 * (qz * qw + qy * qx), sqx - sqy - sqz + sqw);

    euler_ang = [roll, pitch, yaw]';
end
