# Matlab-Delatin
A matlab 3D terrain mesh generation tool. Approximates a height field with a Delaunay triangulation, minimizing the amount of points and triangles for a given maximum error.

Delatin is a port of Michael Fogleman's hmm (C++), which is in turn based on the paper Fast Polygonal Approximation of Terrains and Height Fields (1995) by Michael Garland and Paul Heckbert.

![Raster2MeshGif](https://github.com/AmirGoldental/Matlab-Delatin/blob/main/Raster2Mesh.gif)


Based on [mapbox/delatin](https://github.com/mapbox/delatin) and a priority queue implementation of Andrew Woodward [andrewmwoodward/PriorityQueueMATLAB](https://github.com/andrewmwoodward/PriorityQueueMATLAB).    

Note that I have decided not to use Matlab's functionality for triangulation and developed the edge class instead. The idea is to keep the code generic so that future Julia implementation will be as easy as possible.
