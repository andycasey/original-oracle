

C     FROM KURUCZ ATLAS 9 (OH bf Crossection in cm2)
      REAL*4 FUNCTION OHOP(T,w)
      REAL*8 w,w1
      REAL*4 t
      SAVE w1
      DIMENSION CROSSOH(15,130),CROSSOHT(15)
      DIMENSION C1(150),C2(150),C3(150),C4(150),C5(150)
      DIMENSION C6(150),C7(150),C8(150),C9(150),C10(150)
      DIMENSION C11(150),C12(150),C13(150)
      EQUIVALENCE (CROSSOH(1, 1),C1(1)),(CROSSOH(1,11),C2(1))
      EQUIVALENCE (CROSSOH(1,21),C3(1)),(CROSSOH(1,31),C4(1))
      EQUIVALENCE (CROSSOH(1,41),C5(1)),(CROSSOH(1,51),C6(1))
      EQUIVALENCE (CROSSOH(1,61),C7(1)),(CROSSOH(1,71),C8(1))
      EQUIVALENCE (CROSSOH(1,81),C9(1)),(CROSSOH(1,91),C10(1))
      EQUIVALENCE (CROSSOH(1,101),C11(1))
      EQUIVALENCE (CROSSOH(1,111),C12(1))
      EQUIVALENCE (CROSSOH(1,121),C13(1))
      DATA C1/-30.855,-29.121,-27.976,-27.166,-26.566,-26.106,-25.742,   2.1
     1-25.448,-25.207,-25.006,-24.836,-24.691,-24.566,-24.457,-24.363,   2.1
     2        -30.494,-28.760,-27.615,-26.806,-26.206,-25.745,-25.381,   2.2
     2-25.088,-24.846,-24.645,-24.475,-24.330,-24.205,-24.097,-24.002,   2.2
     3        -30.157,-28.425,-27.280,-26.472,-25.872,-25.411,-25.048,   2.3
     3-24.754,-24.513,-24.312,-24.142,-23.997,-23.872,-23.764,-23.669,   2.3
     4        -29.848,-28.117,-26.974,-26.165,-25.566,-25.105,-24.742,   2.4
     4-24.448,-24.207,-24.006,-23.836,-23.692,-23.567,-23.458,-23.364,   2.4
     5        -29.567,-27.837,-26.693,-25.885,-25.286,-24.826,-24.462,   2.5
     5-24.169,-23.928,-23.727,-23.557,-23.412,-23.287,-23.179,-23.084,   2.5
     6        -29.307,-27.578,-26.436,-25.628,-25.029,-24.569,-24.205,   2.6
     6-23.912,-23.671,-23.470,-23.300,-23.155,-23.031,-22.922,-22.828,   2.6
     7        -29.068,-27.341,-26.199,-25.391,-24.792,-24.332,-23.969,   2.7
     7-23.676,-23.435,-23.234,-23.064,-22.920,-22.795,-22.687,-22.592,   2.7
     8        -28.820,-27.115,-25.978,-25.172,-24.574,-24.115,-23.752,   2.8
     8-23.459,-23.218,-23.017,-22.848,-22.703,-22.579,-22.470,-22.376,   2.8
     9        -28.540,-26.891,-25.768,-24.968,-24.372,-23.914,-23.552,   2.9
     9-23.259,-23.019,-22.818,-22.649,-22.504,-22.380,-22.272,-22.177,   2.9
     A        -28.275,-26.681,-25.574,-24.779,-24.186,-23.729,-23.368,   3.0
     A-23.076,-22.836,-22.636,-22.467,-22.322,-22.198,-22.090,-21.996/   3.0
      DATA C2/-27.993,-26.470,-25.388,-24.602,-24.014,-23.560,-23.200,   3.1
     1-22.909,-22.669,-22.470,-22.301,-22.157,-22.033,-21.925,-21.831,   3.1
     2        -27.698,-26.252,-25.204,-24.433,-23.851,-23.401,-23.043,   3.2
     2-22.754,-22.515,-22.316,-22.148,-22.005,-21.881,-21.773,-21.679,   3.2
     3        -27.398,-26.026,-25.019,-24.267,-23.696,-23.251,-22.896,   3.3
     3-22.609,-22.372,-22.174,-22.007,-21.864,-21.741,-21.634,-21.540,   3.3
     4        -27.100,-25.791,-24.828,-24.102,-23.543,-23.106,-22.756,   3.4
     4-22.472,-22.238,-22.041,-21.875,-21.733,-21.611,-21.504,-21.411,   3.4
     5        -26.807,-25.549,-24.631,-23.933,-23.391,-22.964,-22.621,   3.5
     5-22.341,-22.109,-21.915,-21.751,-21.610,-21.488,-21.383,-21.290,   3.5
     6        -26.531,-25.310,-24.431,-23.761,-23.238,-22.823,-22.488,   3.6
     6-22.214,-21.986,-21.795,-21.633,-21.494,-21.374,-21.269,-21.178,   3.6
     7        -26.239,-25.066,-24.225,-23.585,-23.082,-22.681,-22.356,   3.7
     7-22.089,-21.866,-21.679,-21.520,-21.383,-21.265,-21.162,-21.072,   3.7
     8        -25.945,-24.824,-24.017,-23.405,-22.923,-22.538,-22.223,   3.8
     8-21.964,-21.748,-21.565,-21.410,-21.276,-21.160,-21.059,-20.970,   3.8
     9        -25.663,-24.587,-23.810,-23.222,-22.761,-22.391,-22.088,   3.9
     9-21.838,-21.629,-21.452,-21.300,-21.170,-21.057,-20.958,-20.872,   3.9
     A        -25.372,-24.350,-23.603,-23.038,-22.596,-22.241,-21.950,   4.0
     A-21.710,-21.508,-21.337,-21.190,-21.064,-20.954,-20.858,-20.774/   4.0
      DATA C3/-25.076,-24.111,-23.396,-22.853,-22.429,-22.088,-21.809,   4.1
     1-21.578,-21.384,-21.220,-21.078,-20.957,-20.851,-20.758,-20.676,   4.1
     2        -24.779,-23.870,-23.189,-22.669,-22.261,-21.934,-21.667,   4.2
     2-21.445,-21.259,-21.101,-20.965,-20.848,-20.746,-20.656,-20.578,   4.2
     3        -24.486,-23.629,-22.983,-22.486,-22.095,-21.781,-21.524,   4.3
     3-21.311,-21.132,-20.980,-20.850,-20.737,-20.639,-20.553,-20.478,   4.3
     4        -24.183,-23.382,-22.774,-22.302,-21.928,-21.627,-21.381,   4.4
     4-21.177,-21.005,-20.859,-20.734,-20.625,-20.531,-20.449,-20.376,   4.4
     5        -23.867,-23.127,-22.561,-22.116,-21.761,-21.474,-21.238,   4.5
     5-21.043,-20.878,-20.738,-20.617,-20.513,-20.423,-20.344,-20.274,   4.5
     6        -23.538,-22.862,-22.340,-21.926,-21.592,-21.320,-21.096,   4.6
     6-20.909,-20.751,-20.617,-20.502,-20.402,-20.315,-20.239,-20.172,   4.6
     7        -23.234,-22.604,-22.120,-21.734,-21.422,-21.166,-20.953,   4.7
     7-20.776,-20.625,-20.497,-20.387,-20.291,-20.208,-20.135,-20.071,   4.7
     8        -22.934,-22.347,-21.898,-21.541,-21.250,-21.010,-20.811,   4.8
     8-20.643,-20.500,-20.378,-20.273,-20.182,-20.102,-20.033,-19.971,   4.8
     9        -22.637,-22.092,-21.676,-21.345,-21.075,-20.853,-20.666,   4.9
     9-20.508,-20.374,-20.259,-20.159,-20.073,-19.997,-19.931,-19.872,   4.9
     A        -22.337,-21.835,-21.452,-21.147,-20.899,-20.693,-20.520,   5.0
     A-20.373,-20.247,-20.139,-20.046,-19.964,-19.892,-19.830,-19.774/   5.0
      DATA C4/-22.049,-21.584,-21.230,-20.950,-20.721,-20.531,-20.372,   5.1
     1-20.236,-20.119,-20.019,-19.931,-19.855,-19.788,-19.729,-19.676,   5.1
     2        -21.768,-21.337,-21.011,-20.754,-20.544,-20.370,-20.223,   5.2
     2-20.098,-19.991,-19.898,-19.817,-19.746,-19.683,-19.628,-19.579,   5.2
     3        -21.494,-21.096,-20.796,-20.559,-20.367,-20.208,-20.074,   5.3
     3-19.960,-19.861,-19.776,-19.701,-19.636,-19.578,-19.527,-19.482,   5.3
     4        -21.233,-20.861,-20.585,-20.368,-20.193,-20.048,-19.926,   5.4
     4-19.821,-19.732,-19.654,-19.586,-19.526,-19.473,-19.426,-19.384,   5.4
     5        -20.983,-20.635,-20.380,-20.181,-20.021,-19.889,-19.778,   5.5
     5-19.683,-19.602,-19.531,-19.469,-19.415,-19.367,-19.324,-19.286,   5.5
     6        -20.743,-20.418,-20.182,-19.999,-19.853,-19.733,-19.633,   5.6
     6-19.547,-19.474,-19.410,-19.354,-19.305,-19.261,-19.223,-19.189,   5.6
     7        -20.515,-20.210,-19.991,-19.824,-19.690,-19.581,-19.490,   5.7
     7-19.413,-19.347,-19.290,-19.240,-19.196,-19.157,-19.122,-19.092,   5.7
     8        -20.297,-20.011,-19.808,-19.654,-19.532,-19.434,-19.352,   5.8
     8-19.282,-19.223,-19.172,-19.127,-19.088,-19.054,-19.023,-18.996,   5.8
     9        -20.090,-19.822,-19.633,-19.491,-19.381,-19.291,-19.218,   5.9
     9-19.156,-19.103,-19.057,-19.018,-18.983,-18.952,-18.925,-18.901,   5.9
     A        -19.893,-19.642,-19.467,-19.337,-19.236,-19.155,-19.089,   6.0
     A-19.034,-18.987,-18.946,-18.912,-18.881,-18.854,-18.831,-18.810/   6.0
      DATA C5/-19.705,-19.472,-19.309,-19.190,-19.098,-19.025,-18.966,   6.1
     1-18.917,-18.876,-18.840,-18.810,-18.783,-18.760,-18.739,-18.721,   6.1
     2        -19.527,-19.310,-19.161,-19.051,-18.968,-18.903,-18.851,   6.2
     2-18.807,-18.771,-18.740,-18.713,-18.690,-18.670,-18.653,-18.637,   6.2
     3        -19.357,-19.159,-19.022,-18.922,-18.847,-18.789,-18.743,   6.3
     3-18.704,-18.673,-18.646,-18.623,-18.603,-18.586,-18.571,-18.558,   6.3
     4        -19.195,-19.016,-18.892,-18.803,-18.736,-18.684,-18.643,   6.4
     4-18.610,-18.583,-18.560,-18.540,-18.523,-18.509,-18.496,-18.485,   6.4
     5        -19.042,-18.883,-18.772,-18.693,-18.634,-18.589,-18.553,   6.5
     5-18.525,-18.501,-18.481,-18.465,-18.451,-18.438,-18.428,-18.419,   6.5
     6        -18.894,-18.758,-18.662,-18.593,-18.542,-18.503,-18.473,   6.6
     6-18.448,-18.428,-18.412,-18.398,-18.386,-18.376,-18.367,-18.359,   6.6
     7        -18.752,-18.639,-18.559,-18.501,-18.458,-18.426,-18.400,   6.7
     7-18.380,-18.363,-18.350,-18.338,-18.328,-18.320,-18.313,-18.306,   6.7
     8        -18.611,-18.523,-18.460,-18.415,-18.381,-18.355,-18.334,   6.8
     8-18.318,-18.304,-18.293,-18.284,-18.276,-18.269,-18.263,-18.258,   6.8
     9        -18.471,-18.408,-18.362,-18.329,-18.304,-18.285,-18.269,   6.9
     9-18.257,-18.247,-18.238,-18.231,-18.224,-18.219,-18.214,-18.210,   6.9
     A        -18.330,-18.290,-18.261,-18.239,-18.223,-18.211,-18.201,   7.0
     A-18.192,-18.185,-18.179,-18.174,-18.169,-18.165,-18.162,-18.159/   7.0
      DATA C6/-18.190,-18.168,-18.154,-18.143,-18.135,-18.129,-18.124,   7.1
     1-18.120,-18.116,-18.112,-18.109,-18.106,-18.104,-18.102,-18.100,   7.1
     2        -18.055,-18.047,-18.043,-18.042,-18.040,-18.039,-18.039,   7.2
     2-18.038,-18.037,-18.036,-18.035,-18.034,-18.033,-18.033,-18.032,   7.2
     3        -17.929,-17.931,-17.935,-17.939,-17.943,-17.946,-17.948,   7.3
     3-17.950,-17.952,-17.953,-17.955,-17.956,-17.957,-17.958,-17.959,   7.3
     4        -17.818,-17.826,-17.834,-17.842,-17.849,-17.855,-17.860,   7.4
     4-17.865,-17.869,-17.872,-17.875,-17.878,-17.881,-17.883,-17.886,   7.4
     5        -17.724,-17.736,-17.747,-17.758,-17.767,-17.775,-17.782,   7.5
     5-17.788,-17.793,-17.798,-17.803,-17.807,-17.811,-17.815,-17.819,   7.5
     6        -17.651,-17.665,-17.678,-17.690,-17.701,-17.710,-17.718,   7.6
     6-17.725,-17.732,-17.738,-17.744,-17.749,-17.755,-17.760,-17.765,   7.6
     7        -17.601,-17.615,-17.629,-17.642,-17.653,-17.663,-17.672,   7.7
     7-17.680,-17.688,-17.695,-17.701,-17.708,-17.714,-17.720,-17.726,   7.7
     8        -17.572,-17.587,-17.602,-17.614,-17.626,-17.636,-17.645,   7.8
     8-17.654,-17.662,-17.670,-17.677,-17.684,-17.691,-17.698,-17.704,   7.8
     9        -17.565,-17.581,-17.595,-17.607,-17.619,-17.629,-17.638,   7.9
     9-17.647,-17.656,-17.664,-17.671,-17.679,-17.686,-17.693,-17.700,   7.9
     A        -17.580,-17.594,-17.608,-17.620,-17.630,-17.640,-17.650,   8.0
     A-17.658,-17.667,-17.675,-17.682,-17.690,-17.697,-17.704,-17.711/   8.0
      DATA C7/-17.613,-17.626,-17.639,-17.649,-17.659,-17.669,-17.677,   8.1
     1-17.686,-17.694,-17.701,-17.709,-17.716,-17.723,-17.730,-17.737,   8.1
     2        -17.663,-17.675,-17.685,-17.695,-17.703,-17.711,-17.719,   8.2
     2-17.727,-17.734,-17.741,-17.748,-17.755,-17.761,-17.768,-17.774,   8.2
     3        -17.728,-17.737,-17.745,-17.752,-17.759,-17.766,-17.772,   8.3
     3-17.778,-17.785,-17.791,-17.797,-17.803,-17.808,-17.814,-17.820,   8.3
     4        -17.803,-17.809,-17.814,-17.818,-17.823,-17.828,-17.832,   8.4
     4-17.837,-17.842,-17.847,-17.852,-17.856,-17.861,-17.866,-17.871,   8.4
     5        -17.884,-17.886,-17.888,-17.889,-17.891,-17.893,-17.896,   8.5
     5-17.899,-17.902,-17.905,-17.908,-17.912,-17.915,-17.919,-17.922,   8.5
     6        -17.966,-17.964,-17.961,-17.959,-17.958,-17.958,-17.958,   8.6
     6-17.959,-17.960,-17.961,-17.963,-17.964,-17.966,-17.968,-17.970,   8.6
     7        -18.040,-18.034,-18.028,-18.023,-18.019,-18.016,-18.013,   8.7
     7-18.012,-18.010,-18.010,-18.009,-18.009,-18.009,-18.009,-18.010,   8.7
     8        -18.096,-18.087,-18.078,-18.071,-18.065,-18.059,-18.055,   8.8
     8-18.051,-18.047,-18.045,-18.042,-18.040,-18.039,-18.037,-18.036,   8.8
     9        -18.125,-18.115,-18.105,-18.097,-18.089,-18.082,-18.076,   8.9
     9-18.070,-18.065,-18.061,-18.057,-18.053,-18.051,-18.048,-18.046,   8.9
     A        -18.120,-18.112,-18.103,-18.095,-18.087,-18.079,-18.072,   9.0
     A-18.066,-18.060,-18.055,-18.050,-18.046,-18.042,-18.039,-18.036/   9.0
      DATA C8/-18.083,-18.078,-18.071,-18.064,-18.057,-18.050,-18.044,   9.1
     1-18.037,-18.032,-18.026,-18.022,-18.017,-18.014,-18.010,-18.007,   9.1
     2        -18.025,-18.022,-18.017,-18.012,-18.006,-18.000,-17.994,   9.2
     2-17.989,-17.984,-17.979,-17.975,-17.971,-17.968,-17.965,-17.963,   9.2
     3        -17.957,-17.955,-17.952,-17.948,-17.943,-17.938,-17.934,   9.3
     3-17.929,-17.925,-17.922,-17.918,-17.916,-17.913,-17.911,-17.910,   9.3
     4        -17.890,-17.889,-17.886,-17.882,-17.879,-17.875,-17.871,   9.4
     4-17.867,-17.864,-17.862,-17.860,-17.858,-17.857,-17.856,-17.855,   9.4
     5        -17.831,-17.829,-17.826,-17.822,-17.819,-17.815,-17.812,   9.5
     5-17.810,-17.807,-17.806,-17.804,-17.803,-17.803,-17.803,-17.803,   9.5
     6        -17.786,-17.782,-17.777,-17.773,-17.769,-17.766,-17.763,   9.6
     6-17.761,-17.759,-17.758,-17.757,-17.757,-17.757,-17.758,-17.759,   9.6
     7        -17.753,-17.747,-17.741,-17.735,-17.731,-17.727,-17.724,   9.7
     7-17.722,-17.721,-17.720,-17.720,-17.720,-17.721,-17.722,-17.724,   9.7
     8        -17.733,-17.724,-17.716,-17.709,-17.703,-17.699,-17.696,   9.8
     8-17.694,-17.693,-17.692,-17.692,-17.693,-17.694,-17.695,-17.697,   9.8
     9        -17.723,-17.711,-17.700,-17.691,-17.685,-17.680,-17.676,   9.9
     9-17.674,-17.673,-17.672,-17.673,-17.673,-17.675,-17.676,-17.678,   9.9
     A        -17.718,-17.702,-17.689,-17.679,-17.672,-17.667,-17.663,  10.0
     A-17.660,-17.659,-17.659,-17.659,-17.660,-17.661,-17.663,-17.665/  10.0
      DATA C9/-17.713,-17.695,-17.681,-17.670,-17.662,-17.656,-17.653,  10.1
     1-17.650,-17.649,-17.649,-17.649,-17.650,-17.651,-17.653,-17.655,  10.1
     2        -17.705,-17.686,-17.671,-17.660,-17.652,-17.647,-17.643,  10.2
     2-17.641,-17.640,-17.640,-17.640,-17.641,-17.643,-17.645,-17.647,  10.2
     3        -17.690,-17.671,-17.657,-17.647,-17.640,-17.635,-17.632,  10.3
     3-17.630,-17.630,-17.630,-17.631,-17.632,-17.634,-17.636,-17.639,  10.3
     4        -17.667,-17.649,-17.637,-17.629,-17.623,-17.619,-17.618,  10.4
     4-17.617,-17.617,-17.618,-17.619,-17.621,-17.623,-17.626,-17.628,  10.4
     5        -17.635,-17.621,-17.611,-17.605,-17.601,-17.600,-17.599,  10.5
     5-17.599,-17.601,-17.602,-17.604,-17.607,-17.609,-17.612,-17.615,  10.5
     6        -17.596,-17.585,-17.579,-17.576,-17.575,-17.575,-17.576,  10.6
     6-17.578,-17.580,-17.582,-17.585,-17.588,-17.591,-17.595,-17.598,  10.6
     7        -17.550,-17.544,-17.542,-17.542,-17.544,-17.546,-17.548,  10.7
     7-17.552,-17.555,-17.558,-17.562,-17.566,-17.570,-17.573,-17.577,  10.7
     8        -17.501,-17.500,-17.501,-17.504,-17.508,-17.513,-17.517,  10.8
     8-17.521,-17.526,-17.530,-17.535,-17.539,-17.544,-17.548,-17.553,  10.8
     9        -17.449,-17.452,-17.457,-17.463,-17.470,-17.476,-17.482,  10.9
     9-17.488,-17.493,-17.499,-17.504,-17.509,-17.514,-17.519,-17.524,  10.9
     A        -17.396,-17.403,-17.412,-17.420,-17.429,-17.437,-17.444,  11.0
     A-17.451,-17.458,-17.464,-17.470,-17.476,-17.481,-17.487,-17.492/  11.0
      DATAC10/-17.344,-17.355,-17.366,-17.377,-17.387,-17.396,-17.405,  11.1
     1-17.413,-17.420,-17.427,-17.434,-17.440,-17.446,-17.452,-17.458,  11.1
     2        -17.295,-17.307,-17.321,-17.333,-17.345,-17.355,-17.365,  11.2
     2-17.373,-17.382,-17.389,-17.397,-17.404,-17.410,-17.417,-17.423,  11.2
     3        -17.249,-17.264,-17.278,-17.292,-17.304,-17.316,-17.326,  11.3
     3-17.335,-17.344,-17.352,-17.360,-17.368,-17.375,-17.382,-17.389,  11.3
     4        -17.209,-17.225,-17.241,-17.255,-17.268,-17.280,-17.291,  11.4
     4-17.301,-17.310,-17.319,-17.327,-17.335,-17.343,-17.350,-17.357,  11.4
     5        -17.177,-17.194,-17.210,-17.225,-17.239,-17.251,-17.262,  11.5
     5-17.272,-17.282,-17.291,-17.300,-17.308,-17.316,-17.324,-17.331,  11.5
     6        -17.154,-17.172,-17.189,-17.204,-17.218,-17.230,-17.242,  11.6
     6-17.252,-17.262,-17.272,-17.280,-17.289,-17.298,-17.306,-17.314,  11.6
     7        -17.144,-17.162,-17.179,-17.194,-17.208,-17.220,-17.232,  11.7
     7-17.242,-17.253,-17.262,-17.271,-17.280,-17.289,-17.297,-17.306,  11.7
     8        -17.146,-17.164,-17.181,-17.196,-17.210,-17.222,-17.234,  11.8
     8-17.245,-17.255,-17.265,-17.274,-17.283,-17.292,-17.301,-17.309,  11.8
     9        -17.163,-17.180,-17.197,-17.212,-17.225,-17.237,-17.249,  11.9
     9-17.260,-17.270,-17.280,-17.289,-17.298,-17.307,-17.316,-17.325,  11.9
     A        -17.193,-17.211,-17.227,-17.241,-17.254,-17.266,-17.277,  12.0
     A-17.288,-17.298,-17.308,-17.317,-17.327,-17.336,-17.345,-17.353/  12.0
      DATAC11/-17.239,-17.256,-17.271,-17.284,-17.297,-17.309,-17.320,  12.1
     1-17.330,-17.340,-17.350,-17.359,-17.369,-17.378,-17.387,-17.395,  12.1
     2        -17.299,-17.315,-17.329,-17.342,-17.354,-17.365,-17.376,  12.2
     2-17.386,-17.396,-17.405,-17.415,-17.424,-17.433,-17.442,-17.451,  12.2
     3        -17.373,-17.388,-17.402,-17.414,-17.425,-17.436,-17.446,  12.3
     3-17.456,-17.466,-17.475,-17.484,-17.493,-17.502,-17.511,-17.520,  12.3
     4        -17.462,-17.476,-17.489,-17.500,-17.511,-17.521,-17.531,  12.4
     4-17.541,-17.550,-17.559,-17.569,-17.578,-17.587,-17.595,-17.604,  12.4
     5        -17.567,-17.581,-17.592,-17.603,-17.613,-17.623,-17.632,  12.5
     5-17.641,-17.651,-17.660,-17.669,-17.678,-17.686,-17.695,-17.704,  12.5
     6        -17.689,-17.701,-17.712,-17.722,-17.732,-17.741,-17.750,  12.6
     6-17.759,-17.768,-17.777,-17.786,-17.795,-17.803,-17.812,-17.821,  12.6
     7        -17.829,-17.840,-17.851,-17.860,-17.869,-17.878,-17.887,  12.7
     7-17.896,-17.904,-17.913,-17.922,-17.930,-17.939,-17.948,-17.956,  12.7
     8        -17.988,-18.000,-18.010,-18.019,-18.028,-18.036,-18.045,  12.8
     8-18.053,-18.062,-18.070,-18.079,-18.087,-18.096,-18.104,-18.112,  12.8
     9        -18.171,-18.183,-18.192,-18.201,-18.210,-18.218,-18.227,  12.9
     9-18.235,-18.243,-18.252,-18.260,-18.268,-18.277,-18.285,-18.293,  12.9
     A        -18.381,-18.393,-18.403,-18.413,-18.422,-18.430,-18.438,  13.0
     A-18.447,-18.455,-18.463,-18.471,-18.479,-18.487,-18.495,-18.503/  13.0
      DATAC12/-18.625,-18.638,-18.650,-18.660,-18.669,-18.678,-18.687,  13.1
     1-18.695,-18.703,-18.711,-18.719,-18.726,-18.734,-18.742,-18.750,  13.1
     2        -18.912,-18.929,-18.943,-18.955,-18.966,-18.975,-18.984,  13.2
     2-18.993,-19.001,-19.008,-19.016,-19.023,-19.031,-19.038,-19.045,  13.2
     3        -19.260,-19.283,-19.303,-19.320,-19.333,-19.345,-19.355,  13.3
     3-19.364,-19.372,-19.380,-19.387,-19.394,-19.400,-19.407,-19.413,  13.3
     4        -19.704,-19.740,-19.771,-19.796,-19.816,-19.832,-19.845,  13.4
     4-19.855,-19.863,-19.870,-19.876,-19.882,-19.887,-19.892,-19.897,  13.4
     5        -20.339,-20.386,-20.424,-20.454,-20.476,-20.492,-20.502,  13.5
     5-20.509,-20.513,-20.516,-20.518,-20.520,-20.521,-20.523,-20.524,  13.5
     6        -21.052,-21.075,-21.093,-21.105,-21.114,-21.120,-21.123,  13.6
     6-21.125,-21.126,-21.127,-21.128,-21.130,-21.131,-21.133,-21.135,  13.6
     7        -21.174,-21.203,-21.230,-21.255,-21.278,-21.299,-21.320,  13.7
     7-21.339,-21.357,-21.375,-21.392,-21.408,-21.424,-21.439,-21.454,  13.7
     8        -21.285,-21.317,-21.346,-21.372,-21.395,-21.416,-21.435,  13.8
     8-21.452,-21.468,-21.483,-21.497,-21.511,-21.524,-21.536,-21.548,  13.8
     9        -21.396,-21.429,-21.459,-21.486,-21.511,-21.532,-21.551,  13.9
     9-21.569,-21.585,-21.600,-21.614,-21.627,-21.640,-21.652,-21.663,  13.9
     A        -21.516,-21.549,-21.580,-21.609,-21.635,-21.658,-21.678,  14.0
     A-21.696,-21.713,-21.728,-21.742,-21.755,-21.767,-21.779,-21.790/  14.0
      DATAC13/-21.651,-21.681,-21.711,-21.738,-21.763,-21.785,-21.804,  14.1
     1-21.821,-21.837,-21.851,-21.864,-21.876,-21.887,-21.898,-21.908,  14.1
     2        -21.810,-21.831,-21.853,-21.874,-21.893,-21.910,-21.925,  14.2
     2-21.938,-21.950,-21.961,-21.971,-21.980,-21.989,-21.998,-22.006,  14.2
     3        -22.009,-22.016,-22.026,-22.037,-22.048,-22.058,-22.066,  14.3
     3-22.074,-22.081,-22.088,-22.094,-22.099,-22.105,-22.111,-22.117,  14.3
     4        -22.353,-22.317,-22.296,-22.284,-22.276,-22.270,-22.266,  14.4
     4-22.262,-22.260,-22.258,-22.257,-22.257,-22.257,-22.258,-22.259,  14.4
     5        -22.705,-22.609,-22.552,-22.515,-22.488,-22.468,-22.451,  14.5
     5-22.438,-22.427,-22.418,-22.410,-22.405,-22.400,-22.397,-22.395,  14.5
     6        -22.889,-22.791,-22.731,-22.690,-22.659,-22.634,-22.612,  14.6
     6-22.594,-22.579,-22.566,-22.555,-22.546,-22.539,-22.533,-22.528,  14.6
     7        -23.211,-23.109,-23.041,-22.989,-22.945,-22.906,-22.872,  14.7
     7-22.842,-22.816,-22.793,-22.774,-22.757,-22.743,-22.732,-22.722,  14.7
     8        -25.312,-24.669,-24.250,-23.959,-23.746,-23.587,-23.463,  14.8
     8-23.366,-23.288,-23.225,-23.173,-23.131,-23.095,-23.066,-23.041,  14.8
     9        -25.394,-24.752,-24.333,-24.041,-23.829,-23.669,-23.546,  14.9
     9-23.449,-23.371,-23.308,-23.256,-23.214,-23.178,-23.149,-23.124,  14.9
     A        -25.430,-24.787,-24.369,-24.077,-23.865,-23.705,-23.582,  15.0
     A-23.484,-23.407,-23.344,-23.292,-23.249,-23.214,-23.185,-23.160/  15.0
      DATA w1/0./
   10 OHOP=0.
      IF(w.EQ.w1) GO TO 30
      w1 =w
      WAVENO=1./w
      EVOLT =WAVENO/8065.479
      N =EVOLT*10.-20.
      EN=FLOAT(N)*.1+2.
      IF(N.LE.0)  RETURN
      IF(N.GE.130)RETURN
      DO 21 IT=1,15
   21  CROSSOHT(IT)=(CROSSOH(IT,N)+(CROSSOH(IT,N+1)-CROSSOH(IT,N))*
     * (EVOLT-EN)/.1)
   30 IF(T.GE.9000.)RETURN
      IF(N.LE.0)RETURN
      IF(N.GE.130)RETURN
      IT=(T-1000.)/200.+1.
      IT=MAX0(IT,1)
      TN=FLOAT(IT)*200.+800.
      IT=(T-2000.)/500.+1.
      IT=MAX0(IT,1)
      TN=FLOAT(IT)*500.+1500.
      OHOP=EXP ( ( CROSSOHT(IT) + 
     *            (CROSSOHT(IT+1)-CROSSOHT(IT)) * (T-TN)/500.)
     *             * 2.30258509299405E0)
      RETURN
      END
