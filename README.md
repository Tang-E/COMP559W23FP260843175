# COMP559W23FP260843175
Edwin Pan's COMP559 Fundamentals of Computer Animation Winter 2023 Final Project repository. A Visual Studio C++ adaptation of Matthias Muller's FLIP PIC fluid simulator (found here https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/18-flip.html) with a small (and not very functional) tweak attempting to add cohesion forces.


# There is at least one bug!!!

std::vector causes a crash when attempting to run the not-fast cohesion force calculation for more than a few seconds. I had wanted to make a video comparison of how the fluid would behave different with and without cohesion forces, but unfortunately I was not able to do this because it crashes too quickly to reach a point of comparison.


# Code and Platform Compatibility

Hi! Just a quick heads up that I wrote this code in Visual Studio on my Windows 10 machine and the libararies were included purely through Visual Studio's system. 

I understand that using CMake and such would have made the code more compatible (as was done in the assignments), but I haven't a darn clue how to make that work and was not able to reverse engineer what the prof and the TA's did for the assignments so I ended up doing it according to how the tutorials I used to learn OpenGL (bless be The Cherno) implemented library importing. If this causes issues, I'm really sorry about that. I guess I'm still stuck in the Microsoft ecosystem LUL. Does Visual Studio even run on Linux?

Anyway if there's one thing this assignment has taught me (aside from how cool fluids are!), it's that I'm extremely grateful for tools like Unity, Python, etc. for how easy they make some things happen. Would have certainly made the prospect of trying to get the code to run in 3d a lot easier since I wouldn't need to learn how to do 3D OpenGL from scratch.


# Future Work (as if I'm coming back to this project lmao unless if I get paid)

I only thought about this as I was producing the writeup for this project, but a lot of the cohesion force calculations are not necessary: we run over particles that are completely surrounded by other particles ("deep water", you could say) when we only really need to run the cohesion forces on particles exposed to air cells. A filter would therefore be to skip cells that do not have neighbouring air cells.
