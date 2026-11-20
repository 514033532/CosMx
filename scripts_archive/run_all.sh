#!/usr/bin/env bash
# 批量处理全部切片:04 分析(CSV)→ 05 图像导出
set -euo pipefail
export PANEL2_ROOT="/home/zhh2/YUBO-202508/ranse/panel2_analysis"

python3 - <<'PYEOF'
import json, subprocess
root = "/home/zhh2/YUBO-202508/ranse/panel2_analysis"
qp = "/home/zhh2/software/QuPath-v0.7.0-Linux/QuPath/bin/QuPath"
d = json.load(open(f"{root}/scripts/slides.json"))
for sid, info in d.items():
    print(f"===== 04 {sid} =====", flush=True)
    rc1 = subprocess.run(
        [qp, "script", "-i", info["path"], "-a", f"slide_id={sid}",
         "-a", f"slides_json={root}/scripts/slides.json",
         f"{root}/scripts/04_process_slide.groovy"],
        stdout=open(f"{root}/images/qc/proc_{sid}_stdout.log", "w"),
        stderr=subprocess.STDOUT).returncode
    print(f"----- 04 exit={rc1} -----", flush=True)
    rc2 = subprocess.run(
        [qp, "script", "-i", info["path"], "-a", f"slide_id={sid}",
         f"{root}/scripts/05_export_images.groovy"],
        stdout=open(f"{root}/images/qc/images_{sid}_stdout.log", "w"),
        stderr=subprocess.STDOUT).returncode
    print(f"===== {sid} DONE (04={rc1}, 05={rc2}) =====", flush=True)
PYEOF
