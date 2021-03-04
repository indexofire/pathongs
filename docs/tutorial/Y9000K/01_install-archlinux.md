# 安装 ArchLinux

制作 usb 安装介质，

```bash
$ 
```

将U盘插入，开机按F2进入BIOS，修改启动选项为新制作的U盘安装介质，保存退出BIOS。重启后电脑引导进入 ArchLinux 安装界面。

```bash
# 先激活无线联网
root@arch-iso # wpa_suppliant -i wlan0 -c <(wpa_passphrase "SSID" "PWD") -B

#
``` 

