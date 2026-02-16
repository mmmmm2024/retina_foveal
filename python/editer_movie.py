import subprocess

input_file = "membrane_ON_CBC.mp4"
output_file = "ONCB_test.mp4"
start_time = 3  # 秒

cmd = [
    "ffmpeg",
    "-ss", str(start_time),
    "-i", input_file,
    "-c", "copy",
    output_file
]

subprocess.run(cmd)
