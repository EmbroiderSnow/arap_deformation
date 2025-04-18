# ARAP-Deformation with Visualization

## 本地部署

请使用Ubuntu22.04，项目依赖OpenMesh、Eigen、OpenGL

OpenMesh和Eigen已经包含在项目文件夹中，安装OpenGL：

```sh
 sudo apt-get install libglfw3-dev libglew-dev
 # 确保当前在arap-deformation目录下
 make build
 # 进行测试
 make test
```

## 基本操作

程序运行后，按下S进入选择模式，可以左键选择约束点，再次点击同一个顶点来取消选择

按下D进入变形模式，鼠标拖动-松开后将控制约束点变形的方向（当前变形模式有bug，在变形之后网格会消失）

解决方法：按下E可以在build目录下保存变形后网格为output.obj，执行
```sh
cd build
./arap-deformation ./output.obj
```
来查看变形效果

按下Q退出选择模式或变形模式

按住鼠标右键向上或向下拖动来进行视角缩放，按住鼠标左键可以转动视角

按C清除所有选中的约束点

（ps：其实使用体验很糟糕，但是工程量实在有点大+已经绞尽脑汁了，还希望多担待QwQ）