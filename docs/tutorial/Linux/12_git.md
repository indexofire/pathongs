# Git

{{ git_page_authors }} 更新于: {{ git_revision_date }}

---

!!! Abstract "内容简介"
    本教程简介git命令的使用，以及终端下用户界面lazygit的简介。

## 安装

在 Archlinux 中安装 git

```bash
$ sudo pacman -S git
```

## 配置

```bash
# 设置颜色 
$ git config --global color.ui true

# 设置命令缩写
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

## 使用

基本使用方法

```bash
# 克隆代码仓库
$ git clone https://github.com/linux/linux.git

# 修改代码后提交所有修改到 master 仓库中，新版本 github 默认仓库名为 main
$ git add -A
$ git commit -m 'fix code'
$ git push origin master
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

## Lazygit

不同操作系统有各种git的gui应用程序，个人通常使用lazygit命令来实现git功能。

```bash
# aur 方式安装
$ git clone https://aur.archlinux.org/lazygit
$ cd lazygit
$ makepkg -si

# lg作为运行命令
$ echo "alias lg='lazygit'" >> $HOME/.bashrc
$ source $HOME/.bashrc
```

使用 lazygit

```bash
# 进入代码仓库根目录
$ cd path-to-repos

# 进入lazygit界面
$ lg
```

lazygit用户界面分为左侧5个功能窗口和右侧一个内容窗口。当在功能窗口切换时，右侧内容窗口也会显示当前功能窗口所包含的内容。按数字1-5切换各个功能窗口。

**功能窗口**

1. Status
2. Files
3. Local Branches
4. Commits
5. Stash

**快捷命令

`p`: git pull
`P`: git push
