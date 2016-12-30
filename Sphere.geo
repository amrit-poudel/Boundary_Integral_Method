
// Mesh.ElementOrder = 2 ;

// Mesh.SecondOrderLinear = 1 ; 

Mesh.Algorithm3D = 4 ; // 3D mesh algorithm (1=Delaunay(+Tetgen), 4=Frontal)

// Mesh.Optimize = 1 ; // Optimize the mesh to improve the quality of tetrahedral elements ( for Tetgen)

Mesh.OptimizeNetgen = 1; // Optimize the mesh using Netgen to improve the quality of tetrahedral elements

//Mesh.Algorithm = 6;

// Sphere parameters
// Center
xc  = 0.0;
yc = 0.0;
zc = 0.0;

// Radius
r = 10.0;


// Mesh parameter
lc = 0.6;

// Geometry
Point(1) = {xc, yc, zc, lc};
Point(2) = {xc+r, yc, zc, lc};
Point(3) = {xc, yc+r, zc, lc};
Point(4) = {xc-r, yc, zc, lc};
Point(5) = {xc, yc-r, zc, lc};
Point(6) = {xc, yc, zc-r, lc};
Point(7) = {xc, yc, zc+r, lc};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {3, 1, 6};
Circle(6) = {6, 1, 5};
Circle(7) = {5, 1, 7};
Circle(8) = {7, 1, 3};
Circle(9) = {2, 1, 7};
Circle(10) = {7, 1, 4};
Circle(11) = {4, 1, 6};
Circle(12) = {6, 1, 2};

Line Loop(13) = {2, 8, -10};
Ruled Surface(14) = {13};
Line Loop(15) = {10, 3, 7};
Ruled Surface(16) = {15};
Line Loop(17) = {-8, -9, 1};
Ruled Surface(18) = {17};
Line Loop(19) = {-11, -2, 5};
Ruled Surface(20) = {19};
Line Loop(21) = {-5, -12, -1};
Ruled Surface(22) = {21};
Line Loop(23) = {-3, 11, 6};
Ruled Surface(24) = {23};
Line Loop(25) = {-7, 4, 9};
Ruled Surface(26) = {25};
Line Loop(27) = {-4, 12, -6};
Ruled Surface(28) = {27};
Surface Loop(29) = {28, 26, 16, 14, 20, 24, 22, 18};
Volume(30) = {29};


// Physical surface and volume
Physical Surface(1) = {28, 26, 16, 14, 20, 24, 22, 18};
Physical Volume(2) = {30};


// Color mesh based on physical entity
Mesh.ColorCarousel = 2; // (0=by element type, 1=by elementary entity, 2=by physical entity, 3=by partition)


/*
If physical line surface or volume is not specificed then gmsh will save all 1D, 2D and 3D mesh elements. 
*/


// scuff-em only reads elements with three nodes. Therefore, no need to specify "physical" keyword. 




