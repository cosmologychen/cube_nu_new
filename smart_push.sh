#!/bin/bash

# 确保在 Git 仓库中工作
if ! git rev-parse --is-inside-work-tree > /dev/null 2>&1; then
  echo "错误：当前不在 Git 仓库中！"
  exit 1
fi

# 1. 执行清理操作
echo "🚀 执行 make clean..."
if [ -f "Makefile" ]; then
  make clean && cd utilities/ && make clean && cd .. && rm *.log
  echo "✅ 清理完成"
else
  echo "⚠️  未找到 Makefile，跳过清理"
fi

# 2. 添加所有更改
echo "📦 添加所有更改到暂存区..."
git add -A

# 3. 显示变更状态
echo "📊 当前变更状态："
git status --short

# 4. 获取提交信息
read -p "✏️  请输入提交说明: " commit_msg
if [ -z "$commit_msg" ]; then
  commit_msg="自动提交于 $(date '+%Y-%m-%d %H:%M')"
fi

# 5. 创建提交
echo "💾 创建提交..."
git commit -m "$commit_msg"

# 6. 获取当前分支
current_branch=$(git branch --show-current)
echo "🌿 当前分支: $current_branch"

# 7. 推送到远程
echo "☁️  推送到远程仓库..."
git push origin "$current_branch"

echo "✅ 操作完成！已推送 $current_branch 分支"
