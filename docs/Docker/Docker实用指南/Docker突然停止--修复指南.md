---
title: Docker突然停止--修复指南
description: 非重启的方法解决Docker停止的问题
tags: [command-line]
sidebar_position: 1
---

# Docker突然停止--修复指南

服务器的 Docker 在一天晚上突然停止运行，目前具体原因尚不明确，但修复过程中遇到的问题和解决方法值得记录。主要遇到以下两个问题：
- **Docker 无法重启** 
- **CUDA 掉驱动**

## 1、Docker 无法启动

### 1.1 问题描述
正常情况下，使用 `sudo systemctl restart docker` 命令来重启 Docker 服务。之后使用 `docker ps` 检查 Docker 是否重启完成时遇到问题，导致 Docker 无法重新启动，具体错误信息如下：
```plaintext
Cannot connect to the Docker daemon at unix:///var/run/docker.sock. Is the docker daemon running?
```
经过全网搜索，发现问题出在 `dockerd` 文件上，但该文件被修改的原因目前还不清楚，留待后续研究。
---

### 1.2 解决方案
在我的情况中，Docker 守护进程是通过 `systemctl` 作为服务启动的。可以使用 `systemctl` 命令找到服务文件路径，示例如下（标记 `^^^` 是我添加的，用于指向一行中的位置，并非 shell 输出的一部分）：
```plaintext
sudo systemctl status docker
● docker.service - Docker Application Container Engine
   Loaded: loaded (/usr/lib/systemd/system/docker.service; disabled; vendor preset: disabled)
                     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   Active: active (running) since Fri 2025-02-07 06:17:18 CST; 16h ago
     Docs: https://docs.docker.com
 Main PID: 28733 (dockerd)
    Tasks: 160
   Memory: 161.5G
   CGroup: /system.slice/docker.service
           ├─28733 /usr/bin/dockerd -H unix:// --containerd=/run/containerd/containerd.sock
           └─61323 /usr/bin/docker-proxy -proto tcp -host-ip 0.0.0.0 -host-port 8080 -container-ip 172.17.0.4 -container-port 8080
```
找到 `docker.service` 文件的地址后，切换到该目录并备份文件：
```plaintext
cd /usr/lib/systemd/system/
# 复制一份文件留作备用
cp docker.service docker.service.cp
```
接下来修改 `docker.service` 文件，先查看包含 `ExecStart=/usr/bin/dockerd` 的行：
```plaintext
# 查看 docker.service 中包含“ExecStart=/usr/bin/dockerd”的一行语句
cat /usr/lib/systemd/system/docker.service | grep ExecStart=/usr/bin/dockerd
# 这是结果
ExecStart=/usr/bin/dockerd -H fd:// --containerd=/run/containerd/containerd.sock
```
将 “dockerd” 命令中的 `-H` 参数修改为使用 Unix 套接字，而不是 `fd`。需要修改的地方使用 `^^` 标注：
```plaintext
ExecStart=/usr/bin/dockerd -H unix:///var/run/docker.sock --containerd=/run/containerd/containerd.sock
                           ^^^
```
保存更改后，重新加载守护进程配置并重启 Docker 服务：
```plaintext
sudo systemctl daemon-reload
sudo systemctl restart docker
```
使用 `docker ps` 检查，发现 Docker 已成功启动。
---

## 2、CUDA 掉驱动

### 2.1 问题描述
尝试重启调用 CUDA 的镜像时，发现无法启动，但仅使用 CPU 的镜像可以正常使用，具体报错信息如下：
```plaintext
docker start zhb_acrs_pst
Error response from daemon: failed to create task for container: failed to create shim task: OCI runtime create failed: runc create failed: unable to start container process: error during container init: error running hook #0: error running hook: exit status 1, stdout: , stderr: Auto-detected mode as 'legacy'
nvidia-container-cli: initialization error: nvml error: driver/library version mismatch: unknown
Error: failed to start containers: zhb_acrs_pst
```
尝试使用 `nvidia-smi` 命令检查，出现以下报错：
```plaintext
nvidia-smi 
Failed to initialize NVML: Driver/library version mismatch
NVML library version: 550.144
```
推测是显卡驱动版本不匹配，但在 Docker 重启之前，显卡驱动可以正常运行。由于服务器上还有其他任务在运行，重新安装驱动需要重启服务器，这种方法不太现实。那么是否有不需要重启服务器的解决方法呢？
---

### 2.2 解决方案
此问题是由于内核模块（kernel mod）的 Nvidia 驱动版本未更新导致的。通常情况下，重启机器可以解决该问题。如果因某些原因无法重启，也可以重新加载内核模块，具体步骤如下：
1. **卸载 Nvidia 内核模块**
2. **重新加载 Nvidia 内核模块**

具体执行步骤如下：
1. **终止与 NVIDIA 相关的进程**
首先查看哪些进程正在使用 CUDA：
```plaintext
$ lsof -n -w /dev/nvidia*
COMMAND  PID USER   FD   TYPE  DEVICE SIZE/OFF  NODE NAME
X       4041 root   12u   CHR 195,255      0t0 32915 /dev/nvidiactl
X       4041 root   16u   CHR 195,255      0t0 32915 /dev/nvidiactl
X       4041 root   17u   CHR   195,0      0t0 19698 /dev/nvidia0
X       4041 root   18u   CHR   195,0      0t0 19698 /dev/nvidia0
```
然后终止这些进程：
```plaintext
kill 4041
```

--- 
2. **重新加载内核模块**
运行 `nvidia-smi` 命令，该命令会自动重新加载内核模块，稍等片刻：
```plaintext
$ nvidia-smi
Fri Feb  7 22:02:35 2025       
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 550.144.03             Driver Version: 550.144.03     CUDA Version: 12.4     |
|-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA GeForce RTX 3090        Off |   00000000:98:00.0 Off |                  N/A |
| 40%   46C    P0            107W /  350W |       1MiB /  24576MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
                                                                                         
+-----------------------------------------------------------------------------------------+
| Processes:                                                                              |
|  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
|        ID   ID                                                               Usage      |
|=========================================================================================|
|  No running processes found                                                             |
+-----------------------------------------------------------------------------------------+
```
至此，修复成功！
---
## 3、后续问题

### 3.1 仍无法启动含有 GPU 的 Docker 镜像
尝试正常启动镜像时仍然失败，根据提示，可能是缺少 `nvidia-container`。在 root 权限下安装即可：
```plaintext
$ docker restart zhb_acrs_pst 
Error response from daemon: Cannot restart container zhb_acrs_pst: exec: "nvidia-container-runtime-hook": executable file not found in $PATH
```
```plaintext
$ sudo yum install -y nvidia-docker2

···
Installed:
  nvidia-docker2.noarch 0:2.13.0-1                                                                                                                                                                                                                           

Dependency Installed:
  nvidia-container-toolkit.x86_64 0:1.16.1-1                                                                                 nvidia-container-toolkit-base.x86_64 0:1.16.1-1                                                                                

Complete!
```
再次尝试重启镜像，重启成功！
---
## 参考来源
- [https://superuser.com/questions/1741326/how-to-connect-to-docker-daemon-if-unix-var-run-docker-sock-is-not-available](https://superuser.com/questions/1741326/how-to-connect-to-docker-daemon-if-unix-var-run-docker-sock-is-not-available)
- [https://comzyh.com/blog/archives/967/](https://comzyh.com/blog/archives/967/)