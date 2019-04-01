# Linux 安装

!!! note "内容简介"
    本节将之前一篇在小米 15 pro笔记本上安装ArchLinux的教程作为介绍Linux安装的范例。由于许多的Linux发行版都有图形化的安装界面，只要稍微具备一些计算机知识，加上搜索引擎，安装完成难度不大。但是对于计算机、Linux和生物信息命令的学习来说，还是ArchLinux这种命令行安装界面更有帮助。以此安装作为入门读物，希望对Linux和安装有所帮助。

## 小米 pro 15.6 安装 Archlinux 笔记

- 操作系统: ArchLinux
- WM: i3wm
- DM: none

### 选择小米pro

仔细评估了自己的需求，横向比较了dell xps，华为mate，小米pro等几款笔记本。最终在看到小米新出的游戏本并不适合我之后，果断淘宝了小米 pro 15.6 低配版，5000以内可以拿下。个人的需求是一款和儿子共同使用的笔记本。自己主要做数据分析和开发用，同时教儿子编程，还可以玩玩简单的游戏。加上小米pro极简化的设计，非常符合个人的装逼风格，因此果断入手。曾经在新款小米13和15.6中犹豫过，最终在散热+噪音 vs 便携性中妥协了前者。

### 硬件与软件

硬件配置：

- CPU  i5 8250
- 内存  8G
- 硬盘 256G SSD
- 屏幕 15.6 1920x1080 IPS 72% 色域

软件配置：

- 操作系统 ArchLinux
- 窗口管理 i3-gaps
- 文件管理 ranger
- 音乐播放 mpd
- 视频播放 mpv
- 文本编辑 vim
- 图片查看 w3m
- 文档处理 wps/libreoffice
- 图形图像处理 wine+photoship/gimp
- 屏幕录制视频 ffmpeg
- 游戏 steam

吐槽：机器到手基本满意。屏幕左上角有一点漏光，考虑到我大部分是在终端下工作，基本没有影响。至于做工，总体感觉是模仿mac，但是由于成本受限只能做到大概。比如机器边缘的割手感相当强烈，奇怪怎么很少看到评测里有人提到这一点？由于我是直线控，对这一点目前表示可以接受。还有就是很多说说的B/C面缝隙的问题，如果是完美主义者对这一点会觉得很难受。但是日常使用中确实没有什么影响。使用才是第一位的。

BIOS过于简陋了。极简的外形设计点赞，但是内里也过于简单了，难道定位对象是女性等普通用户？起始小米pro配程序员还是挺合适的。

屏幕的优点是IPS 72% 色域的，确实非常亮丽。不过是镜面屏，这一点不太符合我的需求。另外分辨率是1920x1080，属于常规配置。15寸的4K屏幕确实显示非常细腻，有一种矢量般油画的感觉。但是4k的缺点是支持还不是很理想，我常用的一些 Java 软件必须用1.8的JRE，不支持4k显示，用起来会很麻烦。高清屏（1920x1080）虽然颗粒感还是比较明显，但是软件、投影等支持目前还说还是比4K更便利一些。

### 安装ArchLinux

小米的BIOS界面功能真是简单，看了一下拿到手的是0502版BIOS。在windows10里跑了ssd测试。BIOS 必须设置 supervisor password 才能将修改 secure boot 为 disabled。开机按F2进入BIOS，按F12选择引导设备。

开机按F12，选择U盘启动，进入Archlinux界面。插上一个大容量U盘，用来将原硬盘内容进行备份，以备将来需要恢复时使用。查看一下当前磁盘情况：

```bash
root@archiso~ # lsblk -f
NAME		FSTYPE		LABEL			UUID					MOUNTPOINT
loop0		squashfs											/run/archiso/sfs/airrootfs
sda
-sda1		vfat		BANQ			C873-E5BF
sdb			iso9660		ARCH_201803		2018-03-01-15-13-16-00
-sdb1		iso9660		ARCH_201803		2018-03-01-15-13-16-00	/run/archiso/bootmnt
-sdb2		vfat		ARCHISO_EFI		8688-0060
nvme0n1
-nvme0n1p1	vfat		SYSTEM			E479-B82E
-nvme0n1p2
-nvme0n1p3	ntfs		Windows			3ACC7A8CCC7A4265
-nvme0n1p4	ntfs		Recovery image	C2DE7AF4DE7ADFD9
```

当前版本的linux内核已经包含了驱动，不会出现nvme ssd硬盘找不到的情况。小米pro 256G SSD硬盘上有4个分区，分别时系统引导EFI分区，保留分区，win10分区和系统恢复分区。我想保留win10的recovery区，找一个U盘来备份分区到文件。

```bash
root@archiso~ # mkdir -p /backup
root@archiso~ # mount /dev/sda1 /backup
root@archiso~ # dd if=/dev/nvme0n4 of=/backup/mmipro_recovery.bin bs=64k conv=notrunc,noerror,sync status=progress
15680+0 records in
15680+0 records out
10276604480 bytes (1.0 GB, 908 MiB) copied, 0.697636 s, 1.5 GB/s
```

连接wifi网络，可以直接用wifi-menu，也可以用wpa_supplicant命令连接。连接成功后运行dhcpcd，就可以访问互联网络了。

```bash
# 用wifi-menu
root@archiso~ # wifi-menu
# 用wpa_supplicant
root@archiso~ # wpa_supplicant -i wlp2s0 -c <(wpa_passphrase "YourSSID" "YourKey") -B
# 连接网络
root@archiso~ # dhcpcd
root@archiso~ # ping -c 3 www.archlinux.org
```

考虑了很久要不要加密安装，由于是和孩子一起使用的，有加密狗的话可以限制小孩。但是考虑到可能存在的一些意外情况，安装就采用普通模式吧。对小孩的使用限制就从其他方面着手了。

运行 blkdiscard 清除硬盘，SSD上所有数据都会被清楚，请确认准备完成再执行。分区很简单，用fdisk或者parted都可以，喜欢图形化界面的可以用cgdisk。

```bash
# 擦除SSD数据
root@archiso~ # blkdiscard /dev/nvme0n1
# fdisk分区
root@archiso~ # fdisk /dev/nvme0n1

# 创建gpt分区表
Command (m for help): g
# 创建分区
Command (m for help): n
Partition number (1-128, default 1): 1
First sector (2048-500118158, default 2048): 2048
Last sector, +sector or +size{K,M,G,T,P} (2048-500118158, default 500118158): +512M

# parted 分区
root@archiso~ # parted /dev/nvme0n1

```

挂载`/`和`/boot`分区到/mnt

```bash
# 添加挂载boot的路径
root@archiso~ # mkdir /mnt/boot
root@archiso~ # mount /dev/nvme0n1p2 /mnt
root@archiso~ # mount /dev/nvme0n1p1 /mnt/boot
# 安装基础系统
root@archiso~ # pacstrap /mnt base base-devel
# 生成fstab
root@archiso~ # genfstab -U /mnt >> /mnt/etc/fstab
# chroot到新安装的系统中进行操作
root@archiso~ # arch-chroot /mnt /bin/bash
```

```bash
# 安装心系统所必须的无线网络相关软件
[root@archiso /]# pacman -S vim wpa_supplicant networkmanager
```

创建管理员密码，新建用户mark，并建立主用户目录和密码。将用户mark设置为wheel组，获得管理员权限。

```bash
[root@archiso /]# passwd root
[root@archiso /]# useradd -G -m mark
```

以前版本的grub不支持nvme，不确定新版本或者grub-git是否可以完美支持，这里我们就不用grub引导了，直接用systemd-boot来引导系统。

```bash
# 安装efi
[root@archiso /]# bootctl install
[root@archiso /]# vim /boot/efi/entities/arch.conf




[root@archiso /]# vim /boot/efi/load.conf

```

分区时我们没有建立swap分区，因为在SSD中，更喜欢用文件方式使用swap，更为灵活。

```
# mkdir /swap
# fallocate -l 2G /swap/swapfile
# chmod 600 /swap/swapfile
# mkswap /swap/swapfile
# swapon /swap/swapfile
# echo "/swap/swapfile none swap defaults 0 0" >> /etc/fstab
```



### 安装 teminal 工具 urxvt

从剪贴板中粘贴内容到urxvt中，可以用 shift+insert，或者鼠标中键。从urxvt中复制，直接用鼠标选中即可。复制到其他程序时，只要使用鼠标中键即可。

```
# pacman -S xorg xorg-xinit

libgl
# pacman -S rxvt-unicode urxvt-perls
# pacman -S linux-headers
# pacman -S sudo
# EDITOR=vim visudo
%wheel ALL=(ALL) ALL
```

```
# useradd mark -g wheel -m
# passwd mark
New password:
Retype new password:
```

### 字体设置

```
# 安装字体
$ sudo pacman -S ttf-inconsolata wqy-microhei wqy-microhei-light wyq-zenhei wqy-bitmapfont ttf-font-awesome ttf-arphic-ukai ttf-arphic-uming noto-fonts

$ sudo ln -s /etc/fonts/conf.avail/70-no-bitmaps.conf /etc/fonts/conf.d/
$ sudo unlink /etc/fonts/conf.d/10-scal-bitmap-fonts.conf
```

### 设置无线网络

安装 ArchLinux 过程中，我们用过 wifi-menu 工具连接了无线网络。进入安装完毕的系统，也可以以管理员运行 wifi-menu 连接，再运行 dhcpcd，即可上网。如果要启用网络服务，开机自动运行和连接，需要安装 NetworkManager

```bash
##
$ sudo pacman -S networkmanager
$ sudo systemctl enable NetworkManager

## 命令行下建立 wifi 连接
$ nmcli dev wifi connect "your-wifi-ssid" password "wifi-password"
## netoworkmanager
```

### 设置音频输入输出

ArchLinux 的内核模块已经包含 ALSA，没有声音是因为默认设置成静音状态。控制声音，要安装相应的工具。

```bash
## alsa 工具
$ sudo pacman -S alsa-utils alsa-oss

## alsamixer 配置声音，通过ncurse界面。
## 输出Master下方如果显示是MM，则处于mute静音状态，按m键接触
$ alsamixer

## amixer 命令行声音设置，下面命令解除静音
$ amixer sset Master unmute

## 查看所有 mixer
$ amixer scontrols
```

声音配置文件，系统级的在`/etc/asound.conf`，用户级的在`~/.asound.rc

### 屏幕色彩矫正

ICC是1993年由包括Adobe、Apple、Kodak、Microsoft在内的数十家公司发起并成立的，全称是国际色彩联盟（International Color Consortium），简称ICC。这个组织的主要目标，就是要在各个设备、软件之间形成统一的色彩标准，即ICC标准。最终的目标就是让输入设备（如数码相机）、显示设备（如显示器）、输出设备（如打印机）能够自始至终得保证色彩被准确重现。

ICC借助于一个独立的色彩空间即设备色彩特征描述文件连接空间/参考空间（Profile Connection Space，缩写PCS）作为色彩空间转换的中间色彩空间，通过彩色设备色彩空间（RGB或CMYK）和PCS空间之间的联系为该设备建立设备特征描述文件（Profile），从而实现了对色彩的开放性管理，使得色彩传递不依赖于彩色设备。PCS空间是通过XYZ空间或LAB空间来定义的。ICC标准总共规定了7种特征描述文件Profile，其中包括：3种基本设备的特征描述文件Profile，即输入设备Profile文件、显示设备Profile文件和输出设备Profile文件；4个附加的特征描述文件Profile，即设备链接Profile、色彩空间转换Profile、抽象Profile和被命名色Profile

在 Mac 与 Windows 平台上，名称分别被自动加上 .icc 与 .icm 的配置文件扩展名。.icm 与 .icc 扩展名可以互换。

ArchLinux安装完后，默认的屏幕色温偏冷。可以用Win10的icc文件矫正。如果没有备份原有色彩文件，可以用下面方法来加载。相关内容[参见这里](https://www.reddit.com/r/Xiaomi/comments/7y6o7d/anyone_got_a_xiaomi_notebook_pro_screen/)。

```bash
## 可以用 edid-decode 来查看屏幕参数，小米pro的屏幕是京东方BOE
## 生产的。我的这块是2017年18周出品。不关心参数的话可以跳过这一步
$ git clone https://aur.archlinux.org/edid-decode-git.git
$ cd edid-decode
$ makepkg -si
$ edid-decode < /sys/class/drm/card0-eDP-1/edid
Extracted contents:
header:          00 ff ff ff ff ff ff 00
serial number:   09 e5 47 07 00 00 00 00 12 1b
version:         01 04
basic params:    a5 22 13 78 02
chroma info:     1b bb a6 58 55 9d 26 0e 4f 55
established:     00 00 00
standard:        01 01 01 01 01 01 01 01 01 01 01 01 01 01 01 01
descriptor 1:    9c 3b 80 36 71 38 3c 40 30 20 36 00 58 c2 10 00 00 1a
descriptor 2:    fd 2d 80 0e 71 38 28 40 30 20 36 00 58 c2 10 00 00 1a
descriptor 3:    00 00 00 fe 00 42 4f 45 20 43 51 0a 20 20 20 20 20 20
descriptor 4:    00 00 00 fe 00 4e 56 31 35 36 46 48 4d 2d 4e 36 31 0a
extensions:      00
checksum:        30
EDID version: 1.4
Manufacturer: BOE Model 747 Serial Number 0
Made in week 18 of 2017
Digital display
8 bits per primary color channel
DisplayPort interface
Maximum image size: 34 cm x 19 cm
Gamma: 2.20
Supported color formats: RGB 4:4:4
First detailed timing includes the native pixel format and preferred refresh rate
Display x,y Chromaticity:
  Red:   0.6484, 0.3447
  Green: 0.3339, 0.6162
  Blue:  0.1503, 0.0576
  White: 0.3105, 0.3349
Established timings supported:
Standard timings supported:
Detailed mode: Clock 152.600 MHz, 344 mm x 194 mm
               1920 1968 2000 2230 hborder 0
               1080 1083 1089 1140 vborder 0
               +hsync -vsync
               VertFreq: 60 Hz, HorFreq: 68430 Hz
Detailed mode: Clock 117.730 MHz, 344 mm x 194 mm
               1920 1968 2000 2190 hborder 0
               1080 1083 1089 1120 vborder 0
               +hsync -vsync
               VertFreq: 47 Hz, HorFreq: 53757 Hz
ASCII string: BOE CQ
ASCII string: NV156FHM-N61
Checksum: 0x30 (valid)
EDID block does NOT conform to EDID 1.4!
	Missing name descriptor
	Missing monitor ranges
```

没有备份icc文件的话，可以在[这里](https://www.notebookcheck.net/uploads/tx_nbc2/BOE_CQ_______NV156FHM_N61.icm)下载。通过.xinitrc用户级方式在进入桌面系统并加载。

```bash
$ wget https://www.notebookcheck.net/uploads/tx_nbc2/BOE_CQ_______NV156FHM_N61.icm
$ sudo cp BOE*.icm /usr/local/share/color/icc/BOE_CQ_NV156FHM_N61.icm
$ git clone https://aur.archlinux.org/xcalib.git
$ cd xcalib
$ makepkg -si
$ echo "/usr/bin/xcalib -d :0 /usr/local/share/color/icc/BOE_CQ_NV156FHM_N61.icm" >> ~/.xinitrc
```


### 键盘按键设置

安装完成后，会发现键盘功能键没有启用。用`xev | awk -F'[ )]+' '/^KeyPress/ { a[NR+2] } NR in a { printf "%-3s %s\n", $5, $8 }'`测试发现实际上是Fn+F1~12默认设置成了F1～12功能。这在一些笔记本中BIOS可以改。小米笔记本如果要在F1～12和Fn+F1~F12之间切换，只要按下Fn+Esc即可。

```bash
## 查看键盘按键对应功能
$ xmodmap -pke | less
...
```

F1-F12对应按键：

- F1 121
- F2 122
- F3 123
- F4 232
- F5 233
- F6 133
- F7
- F8
- F9
- F10
- F11 107
- F12 118

控制屏幕亮度，xbacklight不起作用，使用light控制

```bash
$ sudo pacman -S light

## eidt in .config/.i3/config
bindsym XF86MonBrightnessUp exec light -A 10
bindsym XF86MonBrightnessDown exec light -U 10

## 如果想查看当前亮度
$ light -G
```

控制声音，这里使用amixer（pactl for pulseAudio）

```bash
## edit in .config/.i3/config
bindsym XF86AudioRaiseVolume exec amixer set Master 10%+
bindsym XF86AudioLowerVolume exec amixer set Master 10%-
bindsym XF86AudioMute exec amixer set Master toggle
```

避免电源按钮误触，进行以下设置。

```bash
$ sudo echo "HandlePowerKey=ignore" >> /etc/systemd/logind.conf

## edit in .config/i3/config
bindsym --release XF86PowerOff mode "$power"
set $power Say Goodbye? [L]ogout | [S]hutdown | [C]ancel
mode "$power" {
	bindsym l exec i3-msg exit
	bindsym s exec systemctl poweroff
    bindsym c mode "default"
	bindsym Return mode "default"
	bindsym Escape mode "default"
}
```


### 设置触摸板

```bash
$ sudo pacman -S xf86-input-libinput
$ sudo vim /etc/X11/xorg.conf.d/20-touchpad.conf
Section "InputClass"
        Identifier "libinput touchpad"
        Driver "libinput"
        MatchIsTouchpad "on"
        MatchDevicePath "/dev/input/event*"
        Option "Tapping" "on"
        Option "ClickMethod" "clickfinger"
        Option "NaturalScrolling" "true"
EndSection
```

### 安装 i3-gaps

```
$ vim .xinitrc
$ exec i3
```

### 安装polybar

```bash
# 安装无线 wifi 支持
$ sudo pacman -S wireless_tools
# 安装 mpd 支持
$ sudo pacman -S libmpdclient mpd
#
$ git clone https://
$ makepkg -si
```

### 安装常用软件

```
# 安装浏览器
$ sudo pacman -S chromium firefox pepper-flash ttf-liberation

# 安装 git
$ sudo pacman -S git-core

# 安装 python3 以及 python 虚拟环境
$ sudo pacman -S python python-pipenv
```



### 安装文件管理器

我们用 ranger 来做为文件管理器。

```
$ pacman -S ranger
```

### 中文输入

```bash
## 安装 fcitx
$ sudo pacman -S fcitx
$ vim ~/.xinitrc
export GTK_IM_MODULE=fcitx
export QT_IM_MODULE=fcitx
export XMODIFIERS=@im=fcitx

## 安装 im 和拼音输入法
$ sudo pacman -S fcitx-im
$ sudo pacman -S fcitx-googlepinyin
```

### 科学上网

修复 sslocal 与 openssl 1.1 冲突，修复方式参见[这里](http://www.54it.top/archives/6224.html)

```bash
$ cd /usr/lib/python3.6/site-packages/shadowsocks/crypto/
$ sudo cp openssl.py openssl.py.bak
$ sudo sed 's/EVP_CIPHER_CTX_cleanup/EVP_CIPHER_CTX_reset/g' openssl.py > openssl.py
```

使用shadowsocks，设置开机自动运行

```bash
## 用户级 shadowsocks
$ sudo sslocal -c ~/.config/shadowsocks/server.json -d start --user mark

## 设置开机自动运行系统级shadowsocks，所有人都可以使用
$ sudo mkdir /etc/shadowsocks
$ sudo ln -s ~/.config/shadowsocks/server.json /etc/shadowsocks/
$ sudo vim /etc/systemd/system/sslocal.service
[Unit]
Description=Daemon to start Shadowsocks Client
Wants=network-online.target
After=network.target
[Service]
Type=simple
ExecStart=/usr/bin/sslocal -c /etc/shadowsocks/server.json --pid-file /var/run/sslocal.pid --log-file /var/log/sslocal.log   
[Install]
WantedBy=multi-user.target

## systemd 设置
$ sudo systemctl enable sslocal.service
$ sudo systemctl start sslocal.service
$ sudo systemctl stop sslocal.service
```

设置浏览器插件

```bash
## SwitchyOmega for Chromium
$ wget https://www.switchyomega.com/static/file/v2.5.6/SwitchyOmega_Chromium.crx
$ 7z x SwitchyOmega_Chromium -r
chromium://extensions
Developer mode
load extensions packs

## SwitchOmega for Firefox
$ wget https://www.switchyomega.com/static/file/v2.5.6/SwitchyOmega_Firefox.fx.xpi
about:addons
Install Add-ons From File
```

本地开启bbr，上传数据（如视频）防阻塞。

```bash
## 开启 bbr
$ sudo echo "net.core.default_qdisc=fq" >> /etc/sysctl.d/233-bbr.conf
$ sudo echo "net.ipv4.tcp_congestion_control=bbr" >> /etc/sysctl.d/233-bbr.conf
$ lsmod | grep bbr
tcp_bbr			20480	1
```

### 双显卡设置与游戏

小米pro15.6双显卡分别是HD620和MX150。后者有2G显存。电脑还是以工作为主，偶尔和儿子玩一玩游戏。小米pro 15.6的BIOS里没有可以选择的显卡，因此无法单独禁用某个显卡。我们采用双显卡用optimus交火来在应用和能耗中寻找平衡。

```bash
$ sudo pacman -S bumblebee
# 添加当前用户到组 bumblebee
$ sudo gpasswd -a `whoami` bumblebee
$ sudo systemctl enable bumblebeed.service

# test bumblebee
$ optirun glxgears -info
```


```bash
## 安装 steam
$ vim /etc/pacman.conf
[multilib]
Include = /etc/pacman.d/mirrorlist

$ sudo pacman -Syu
$ sudo pacman -S steam
```

右键点击游戏-> Properties，添加参数 primusrun


### 电源管理

```bash
$ sudo pacman -S tlp
```

### 设置用户目录路径

```bash
$ sudo pacman -S xdg-user-dirs
$ vim ~/.config/user-dirs.dirs

XDG_DESKTOP_DIR="$HOME/desktop"
XDG_DOCUMENTS_DIR="$HOME/documents"
XDG_DOWNLOAD_DIR="$HOME/downloads"
XDG_MUSIC_DIR="$HOME/musics"
XDG_PICTURES_DIR="$HOME/pictures"
XDG_PUBLICSHARE_DIR="$HOME/public"
XDG_TEMPLATES_DIR="$HOME/templates"
XDG_VIDEOS_DIR="$HOME/videos"

$ xdg-user-dirs-update
```


### 办公软件 wps

```bash
$ sudo echo "[archlinuxcn]" >> /etc/pacman.conf
$ sudo echo "Server = http://repo.archlinuxcn.org/$arch" >> /etc/pacman.conf
$ sudo pacman -Syu
$ sudo pacman -S archlinuxcn-keyring
$ sudo pacman -S wps-office ttf-wps-fonts

# $ git clone https://aur.archlinux.org/ttf-ms-fonts.git
# $ cd ttf-ms-fonts && makepkg -si
```



### screen session

平时工作经常要通过ssh访问服务器，因此使用screen是避免不了的。

### 蓝牙

通过配置，大约在中断蓝牙连接30秒后才能自动链接上。

```bash
$ sudo pacman -S bluez bluez-utils
$ sudo modprobe btusb
$ sudo systemctl enable bluetooth.service
# 添加当前用户到组lp
$ sudo gpasswd -a `whoami` lp

# 启动 bluetoothd
$ sudo systemctl start bluetooth.service

# 确认 hci0 是否被rfkill关闭
$ rfkill list

# 链接 blue 设备
$ bluetoothctl
Agent registered
[bluetooth]

# 开机蓝牙自动开启，并且自动连接之前握手成功的蓝牙鼠标。
$ sudo vim /etc/bluetooth/main.conf
Name=BlueZ
Class=0x0000010c
DiscoverableTimeout=30
Cache=always
ReconnectUUIDs=
RecommectAttempts=7
ReconnectIntervals=1,2,4,8,16,32,64
AutoEnable=true
```




### Conda and Bioconda

```bash
$ cd /tmp
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ sh Miniconda3-latest-Linux-x86
...
# change path
[/home/mark/.conda]

$ conda config --add channels defaults
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
```


### Music 播放

命令行下的音乐播放器众多，比较符合个人喜欢的是cmus/nm

```bash

```



### 指纹识别

用`fprintd`验证时，需要将右手食指划第一节慢慢滑过sensor（2秒左右）5次。

```bash
## Archlinux 源里的版本还不支持 Elan。github下载最新源代码
## 支持了小米pro 04f3:0c1a 的 Elan Microeletronics 指纹识别器。
$ git clone https://github.com/iafilatov/libfprint
$ cd libfprint
$ ./autogen.sh
$ ./autoconfigure
$ make && sudo make install
$ sudo ln -s /usr/local/lib/libfprint.so.0 /usr/lib

## 不安装 libfprint 依赖
$ sudo pacman -Sdd fprintd

## 确认fprintd.service正常运行，如果failed，用`systemctl status fprintd`查看问题
$ systemctl | grep fprinted.service

## 修改登录验证机制，如果顺序是fprintd在前，则登录时先用finger验证
$ sudo vim /etc/pam.d/system-local-login
auth      sufficient pam_fprintd.so
auth      include    system-login

## 建立指纹并验证
$ fprintd-enroll
Using device /net/reactivated/Fprint/Device/0
Enrolling right-index-finger finger.
Enroll result: enroll-stage-passed
Enroll result: enroll-stage-passed
Enroll result: enroll-stage-passed
Enroll result: enroll-stage-passed
Enroll result: enroll-completed

## 验证过程，看到veryify-match即表示验证通过
$ fprintd-verify
Using device /net/reactivated/Fprint/Device/0
Listing enrolled fingers:
 - #0: right-index-finger
Verify result: verify-match (done)

## 重启，由于我没有用 Display manager，用xinit startx来启动i3，
## 所以console界面输入用户名后就会提示“滑动手指验证指纹登录”
$ reboot

## i3lock指纹解锁有一个问题在于需要输错1次密码后才能正常指纹解锁。
## 可能需要 xss-lock 等其他方法解决
$ sudo pacman -S i3lock
$ sudo vim /etc/pam.d/i3lock
auth include login

## 总结：
## 小米的 Elan 指纹识别器面积太小，导致libfprint需要滑动才能识别,
## 另外感觉指纹识别速度不是特别敏锐，扫完指纹需要0.5s左右才能解锁。
## 总体上来说，可用性不是太高。个人感觉不如摄像头人脸识别可用性更高...
```

i3lock 参见[这里](https://www.reddit.com/r/i3wm/comments/7p07nx/not_working_i3lock_with_a_fingerprint_scanner/)

### 摄像头

摄像头默认配置好可以直接使用

```bash
## 用ffmpeg测试
$ ffmpeg -f v4l2 -video_size 640x480 -i /dev/video0 output.mkv
$ mpv output.mkv
```


### Wine

```bash
$ sudo pacman -S wine wine-mono wine-tricks wine_gecko
$ WINEPREFIX=$HOME/.wine wineboot -u
$ WINEARCH=win64 WINEPREFIX=$HOME/.wine64 winetricks fontsmooth=rgb gdiplus vcrun2008 vcrun2010 vcrun2012 vcrun2013 vcrun2015 atmlib msxml3 msxml6 gdiplus corefonts
```

### Photoshop

Adobe Photoshop CC 2018 v19.1.3.49649

```bash
$ wget http://209.126.105.33/DownloadBull/Portable_Adobe_Photoshop_CC_2018_v19.0.0.165_x64.zip
$ 7z e Portable_Adobe_Photoshop_CC_2018_v19.0.0.165_x64.zip
$ wine64 PhotoshopPortable.exe
```


### 美化字体设置

许多发行版如Ubuntu优化过字体渲染效果，而 ArchLinux 从头安装的话需要自己来配置。对于普通人来说这是一个复杂而痛苦的过程，那么最简便的方法就是直接安装 Infinality 套装获得现成的字体优化渲染设置。[Infinality]() 是一个 Freetype 的修改版，在 Linux 下提供更好的字体渲染效果，并拥有简便的定制特性。它设置了其他操作系统的字体风格6，可以方便的切换。

安装 fontconfig-infinality-ultimate

```bash
## 安装 freetype2-infinality
$ sudo pacman -R freetype2
$ gpg --keyserver pgp.mit.edu --recv-keys C1A60EACE707FDA5
$ git clone https://aur.archlinux.org/packages/freetype2-infinality/
$ cd freetype2-infinality && makepkg -si
```

```bash
## 安装 fontconfig-infinality-ultimate
## 先卸载 fontconfig
$ sudo pacman -R fontconfig
$ git clone https://aur.archlinux.org/fontconfig-infinality-ultimate.git
$ cd fontconfig-infinality-ultimate && makepkg -si
```




###屏幕保护和锁屏

如果用 xscreensaver 构建

```bash
bindsym $mod+F12 exec i3lock
```


### 网络工具

ArchLinux 的 net-tools 已经默认 deprecated，原来常用 netstat 命令有[新代替工具](http://inai.de/2008/02/19)。


### vim && Python IDE

```bash
$ mv ~/.vim ~/.vim.old
$ sudo pacman -S python-neovim
$ sh -c "$(curl -fsSL https://raw.githubusercontent.com/jarolrod/vim-python-ide/master/setup.sh)"
$
```

```bash
## 安装 nerd 字体
$ git clone https://aur.archlinux.org/nerd-fonts-complete.git
$ cd nerd-fonts-complete && makepkg -si
$ git clone https://github.com/chriskempson/base16-shell.git ~/.config/base16-shell
$ echo 'BASE16_SHELL=$HOME/.config/base16-shell/' >> ~/.bashrc
$ echo '[ -n "$PS1" ] && [ -s $BASE16_SHELL/profile_helper.sh ] && eval "$($BASE16_SHELL/profile_helper.sh)"' >> ~/.bashrc
```

```
## add these line into .vimrc
if filereadable(expand("~/.vimrc_background"))
  let base16colorspace=256
  source ~/.vimrc_background
endif
```

```bash
## 安装 python git module
$ sudo pacman -S python-gitpython
```

### 图标美化

```bash
$ sudo pacman -S numix-icon-theme-git numix-circle-icon-theme-git numix-gtk-theme
```

### 动态壁纸

Linux 下可以用 xwinwraper 之类的工具将视频或 Gif 动态图片作为壁纸，达到动态桌面的效果。当然代价是牺牲一些CPU和内存。

1. 设置compton
2. 设置xwinwrap
3. 准备视频文件：建议使用延时摄影或

```bash
$ youtube-dl --proxy `socks server` -F `URL`
$ youtube-dl --porxy `socks server` -f 137 `URL`
```


```bash
$ touch run_wallpaper.sh

# run_wallpaper.sh
#!/bin/sh
mpv $1 --wid="$2" --loop
```

```bash
$ xwinwrap -fs -argb -ni -un -s -b -st -sp -nf -ov -- run_wallpaper.sh video.mp4 WID
```

```bash
$ cmatrix
$ ffmpeg -f x11grab -s 1920x1080 -i $DISPLAY output.mp4
```

### 打印机

使用CUPS管理打印机。


```bash
$ lp my.doc
```


1. [Infinality]: http://www.infinality.net/
