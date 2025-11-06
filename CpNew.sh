if [ "$#" -lt 1 ]; then
  echo "No argument provided"
  exit 1
fi

# 检查路径是否为'./xxx'格式
if [[ $1 == './'*'/'* ]]; then
  # 提取路径部分，不包括最后的文件名
  path=$(dirname "$1")
  path=${path#*/}
fi

echo cp -v /home/ChenBH/cube_nu_new/$1 ./$path