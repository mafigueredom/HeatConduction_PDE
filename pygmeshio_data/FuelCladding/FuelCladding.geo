// This code was created by pygmsh vunknown.
SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMin = 0.03;
Mesh.CharacteristicLengthMax = 0.03;
s1 = news;
Disk(s1) = {0.0, 0.0, 0.0, 4.08};
s2 = news;
Disk(s2) = {0.0, 0.0, 0.0, 4.65};
bo1[] = BooleanDifference{ Surface{s2}; Delete; } { Surface{s1}; Delete;};
Physical Curve(20) = {s1};
Physical Curve(30) = {s2};
Physical Surface(200) = {bo1[]};