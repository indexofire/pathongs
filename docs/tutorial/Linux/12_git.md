# Git

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

### 配置

```bash
$ git config --global color.ui true
$ git config --global alias.co checkout
$ git config --global alias.ci commit
$ git config --global alias.st status
$ git config --global alias.br branch
# 设置 Editor 使用 sublime
$ git config --global core.editor "subl -w"
# 设置 Editor 使用 vim
$ git config --global core.editor "vim"
# 列举所有配置
$ git config -l
```

### 比较版本差异

```bash
# 工作区与缓存区的差异
$ git diff
# 缓存区与最新 commit 之间的差异
$ git diff -–cached
$ git diff HEAD # 工作区与最新 commit 之间的差异
$ git diff @ # 同上
```

### 比较当前版本和之前版本，以下三者皆可

```bash
$ git diff @~..@
$ git diff HEAD^ HEAD
$ git show @
$ diff -c file1 file2 # 比较两个文件，上下文模式
$ diff -u file1 file2 # 比较两个文件，统一模式
```

### 撤销修改

```bash
# 撤销工作区编辑
$ git config --global alias.co checkout # git checkout 比较常用，设个别名
$ git co .
$ git co -- filename

# 撤销缓存区
$ git reset HEAD filename # add 之后，恢复后回到编辑状态。
$ git reset @ filename # 同上

# 删除当前历史
$ git log --pretty=oneline # 查看最近的日志
$ git reset --hard @~ # 直接回退到上一个版本
$ git reset --hard HEAD^ # 同上

# 比较安全的做法
$ git revert <$id>    # 恢复某次提交的状态，恢复动作本身也创建了一次提交对象，需重写 commit
$ git revert HEAD     # 恢复最后一次提交的状态

# 感觉自己这次写错了
$ get reset @~ # 回到工作区的编辑状态，跟上方的 hard 还是有区别的。

# 改变git log到标准输出
git config ager.log false
```
