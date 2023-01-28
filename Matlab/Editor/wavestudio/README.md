## wavestudio
A sound editor that's lightweight and performant.

### More Info
Written in modern C++ and uses ImGui to render the user interface. This project is aimed towards providing a 
lightweight program which will be cross platform.

### FAQ
**Why not use Qt Widgets or GTK?**

Qt's API's are proven and effective, but Qt's libraries are often 
bigger than necessary. GTK is difficult to use.

**Why make another audio editor?**

Self improvement. 

### Authors
Ben Cottrell

### Building
This has been tested on Debian 11.
See CMakeLists.txt, I've set a few variables for paths to various components including:
 - ImGui
 - ImGui Addons

Change these variable values to where these components are located.
Make sure SDL2 and an OpenGL library is installed.
