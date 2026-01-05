# source module_load_t7920.sh

echo "编译FoF工具..."
cd utilities/
if ! make fof.x; then
    echo "错误：FoF编译失败！"
    exit 1
fi
cd ..
# 检查编译结果
if [ ! -f "utilities/fof.x" ]; then
    echo "错误：FoF可执行文件未生成！"
    exit 1
fi
echo "FoF编译成功，提交作业..."

qsub run.qsub
# tail -f run.log