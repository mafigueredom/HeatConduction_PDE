// This code was created by pygmsh vunknown.
SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.05;
s0 = news;
Disk(s0) = {0.0, 0.0, 0.0, 4};
Physical Curve(10) = {s0};
Physical Surface(100) = {s0};