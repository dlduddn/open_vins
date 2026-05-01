function yamlName = searchYaml(cfg, t_cam)
    % 목적: t_cam에 대응하는 YAML 파일명을 찾는다.
    % 과정:
    %   1. YAML 파일명은 nanosecond 단위 timestamp라고 가정한다.
    %   2. t_cam은 second 단위 timestamp이므로 millisecond key로 변환한다.
    %   3. YAML 파일명 timestamp도 millisecond key로 변환한 뒤 같은 key가 있는지 비교한다.
    
    persistent cache % YAML 폴더를 매번 스캔하지 않기 위한 함수 내부 cache

    % 최초 호출이거나 검색 대상 YAML 폴더가 바뀐 경우에만 cache를 다시 만든다.
    if isempty(cache) || ~strcmp(cache.yamlDir, cfg.yamlDir)
        cache = buildYamlCache(cfg.yamlDir);
    end

    % camera timestamp를 second -> millisecond key로 변환한다.
    % floor를 사용해 ms보다 작은 자릿수는 버린다.
    camMsKey = int64(floor(t_cam * 1e3));

    % camera timestamp key와 같은 YAML timestamp key를 직접 찾는다.
    % timestamp는 유일해야 하므로 2개 이상 발견되면 데이터/파일명 중복 오류로 처리한다.
    yamlIdxForCam = find(cache.yamlMsKey == camMsKey);
    if numel(yamlIdxForCam) > 1
        error('Duplicate YAML timestamp key matched for t_cam %.9f: msKey=%d', t_cam, camMsKey);
    end
    matchedCam = ~isempty(yamlIdxForCam);

    % 이번 호출에서 매칭된 YAML 파일을 logical vector로 표시한다.
    validYaml = false(cache.nYaml, 1);
    if matchedCam
        validYaml(yamlIdxForCam) = true;
    end

    % 누적 매칭 상태를 갱신한다.
    % 몇 개의 YAML 파일이 매칭되었는지 진행 상황을 확인
    cache.cumulativeValidYaml = cache.cumulativeValidYaml | validYaml;
    nMatched = nnz(cache.cumulativeValidYaml);

    % 매칭된 camera timestamp에 대응하는 YAML 파일명을 반환한다.
    % cache.yamlNames에는 확장자를 제거한 파일명만 저장되어 있다.
    if matchedCam
        yamlName = cache.yamlNames{yamlIdxForCam};
    else
        yamlName = '';
    end

    if getLogical(cfg, 'yamlLogProgress', true)
        fprintf("  yaml  = %s\n", formatNameForDisplay(yamlName));
        fprintf('  %d/%d YAML matched\n\n', nMatched, cache.nYaml);
    end
end

function cache = buildYamlCache(yamlDir)
    % 목적: yamlDir 안의 모든 .yaml 파일을 읽어 timestamp 검색용 cache를 만든다.
    yamlFiles = dir(fullfile(yamlDir, '*.yaml'));
    nYaml = numel(yamlFiles);

    % yamlNames: 확장자를 제거한 YAML 파일명
    % yamlMsKey: 파일명 timestamp를 millisecond 단위로 변환한 key
    yamlNames = cell(nYaml, 1);
    yamlMsKey = zeros(nYaml, 1, 'int64');

    % 각 YAML 파일명에서 timestamp를 읽어 ms key로 변환한다.
    for i = 1:nYaml
        [~, name] = fileparts(yamlFiles(i).name);
        yamlNames{i} = name;
        yamlMsKey(i) = yamlNameToMsKey(name);
    end

    % timestamp 오름차순으로 정렬한다.
    % yamlMsKey와 yamlNames는 같은 sortIdx를 적용해야 key와 파일명이 계속 대응된다.
    [yamlMsKey, sortIdx] = sort(yamlMsKey);

    % searchYaml에서 재사용할 정보를 하나의 구조체로 묶는다.
    cache = struct();
    cache.yamlDir = yamlDir;
    cache.yamlNames = yamlNames(sortIdx);
    cache.yamlMsKey = yamlMsKey;
    cache.nYaml = nYaml;

    % YAML별 누적 매칭 여부를 저장한다.
    % searchYaml이 호출될 때마다 OR 연산으로 갱신된다.
    cache.cumulativeValidYaml = false(nYaml, 1);
end

function msKey = yamlNameToMsKey(name)
    % 목적: nanosecond timestamp 파일명을 millisecond key로 변환한다.

    % [Validate] YAML 파일명은 숫자로만 구성된 nanosecond timestamp여야 한다.
    if ~all(isstrprop(name, 'digit')) || numel(name) < 7
        error('Unexpected YAML filename format: %s', name);
    end

    % ns 단위 전체를 double timestamp와 직접 비교하면 부동소수점 오차로 매칭이 실패할 수 있다.
    % 따라서 마지막 6자리(ns -> ms)를 제거하고, ms 해상도 key만 비교한다.
    msKey = int64(str2double(name(1:end-6)));
end

function formattedName = formatNameForDisplay(Name)
    formattedName = regexprep(char(Name), '(\d{3})$', '($1)');
end
