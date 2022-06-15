
% comp. coeffs
Cxu = -1.*[-1, -0.36, -0.19, -0.11, -0.06, -0.03, 0.00, 0.02, 0.04, 0.05, 0.06, 0.08, 0.08, 0.09, 0.10, 0.11, 0.11, 0.12, 0];
Cxl = -1.*[0, -0.36, -0.19, -0.11, -0.06, -0.03, 0.00, 0.02, 0.04, 0.05, 0.06, 0.08, 0.08, 0.09, 0.10, 0.11, 0.11, 0.12, 0];
Cyu = -1.*[0, 0.93, 0.98, 0.99, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.99, 0.99, 0.99, 0];
Cyl = -1.*[0, -0.932802908, -0.982185387, -0.993704522, -0.998019745, -0.999621206, -0.999999656, -0.999782918, -0.999261187, -0.998609885, -0.997890692, -0.997160784, -0.996452891, -0.995741964, -0.995032303, -0.994319925, -0.9935531, -0.992714591, 0];

% coeffs for chordwise and chordnormal graphs
Cnx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
Cny = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
Cwx = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
Cwy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

% upper curve
Xup = [0, 3.5, 8.8, 14.0, 19.3, 24.6, 29.8, 35.1, 40.4, 45.6, 50.9, 56.2, 61.4, 66.7, 72.0, 77.2, 82.5, 87.8, 100];
Yup = [0, 3.04, 4.46, 5.24, 5.70, 5.93, 6.00, 5.95, 5.79, 5.55, 5.24, 4.87, 4.44, 3.98, 3.47, 2.92, 2.34, 1.72, 0];

%lower cp comp
Uup_0 = [1.511309524, -1.474435196, -1.43010113, -1.420143027, -1.4032162, -1.378346222, -1.345032719, -1.312722949, -1.279429251, -1.247324614, -1.220238095, -1.194163192, -1.167657551, -1.14608076, -1.117122473, -1.022010708, -0.950118765, -0.938058368, 0];
Vup_0c = Uup_0.*Cyu;
Uup_0c = Uup_0.*Cxu;

Uup_4 = [0.725648415, -1.438616715, -1.091276719, -0.932525952, -0.804260219, -0.689516129, -0.621683968, -0.565092166, -0.51384083, -0.458884416, -0.41037464, -0.367570687, -0.323716099, -0.278546713, -0.235870819, -0.193194925, -0.142363112, -0.084246971, 0];
Vup_4 = Uup_4.*Cyu;
Uup_4 = Uup_4.*Cxu;

Uup_6 = [0.29294653, -2.008527572, -1.445769449, -1.007959068, -0.910997732, -0.823094425, -0.74059293, -0.656267725, -0.584659091, -0.513636364, -0.450255827, -0.395904437, -0.33996589, -0.285471056, -0.232254401, -0.181147076, -0.125568182, -0.065946561, 0];
Vup_6 = Uup_6.*Cyu;
Uup_6 = Uup_6.*Cxu;

Uup_8 = [-0.14084507, -2.672655811, -1.304568528, -1.160292463, -1.016901408, -0.903207653, -0.803260259, -0.702641934, -0.618806306, -0.538807649, -0.465953855, -0.403153153, -0.341586944, -0.281761717, -0.224662162, -0.170786517, -0.115233277, -0.061797753, 0];
Vup_8 = Uup_8.*Cyu;
Uup_8 = Uup_8.*Cxu;

Uup_10 = [-0.575027996, -3.168903803, -1.546831183, -1.264541387, -1.086640581, -0.952434247, -0.837898267, -0.723713647, -0.628842929, -0.541596873, -0.463087248, -0.394295302, -0.329234209, -0.267337808, -0.212290503, -0.1640625, -0.119686801, -0.082168809, 0];
Vup_10 = Uup_10.*Cyu;
Uup_10 = Uup_10.*Cxu;

Uup_12 = [0.544134727, -0.718768158, -0.765116279, -0.744606414, -0.752325581, -0.736596737, -0.728596389, -0.751599767, -0.756976744, -0.732132481, -0.703033839, -0.701396973, -0.679464805, -0.655393586, -0.623906706, -0.59895227, -0.565192084, -0.549156486, 0];
Vup_12 = Uup_12.*Cyu;
Uup_12 = Uup_12.*Cxu;

Uup_14 = [0.56982911, -0.611764706, -0.631176471, -0.627058824, -0.619187758, -0.618235294, -0.614479105, -0.610849057, -0.612617925, -0.614976415, -0.617560401, -0.624188791, -0.629345905, -0.638235294, -0.647232038, -0.654684738, -0.665097116, -0.657041839, 0];
Vup_14 = Uup_14.*Cyu;
Uup_14 = Uup_14.*Cxu;

% lower curve
Xlo = Xup;
Ylo = -1*Yup;

% lower cp comp
Ulo_0 = Cxl.*Uup_0;
Vlo_0 = Cyl.*Uup_0;

Ulo_4 = [0.5982958, 0.03406326, -0.135892748, -0.189287888, -0.233151184, -0.256379101, -0.259416768, -0.255319149, -0.242736077, -0.226060606, -0.210686096, -0.195255474, -0.176899696, -0.155528554, -0.134101942, -0.111111111, -0.087537994, -0.069781553, 0];
Vlo_4 = Ulo_4.*Cyl;
Ulo_4 = Ulo_4.*Cxl;

Ulo_6 = [-0.352024922, 0.496577474, 0.211478478, 0.088930348, 0.001246883, -0.053549191, -0.081467662, -0.102053516, -0.110211706, -0.108898569, -0.108343711, -0.106409459, -0.099626401, -0.090342679, -0.08089608, -0.068535826, -0.049129353, -0.02615193, 0];
Vlo_6 = Ulo_6.*Cyl;
Ulo_6 = Ulo_6.*Cxl;

Ulo_8 = [-0.791222571, 0.620451694, 0.316844083, 0.17565872, 0.075862069, 0.011904762, -0.024451411, -0.05339196, -0.067584481, -0.072818581, -0.077791719, -0.080877743, -0.079524108, -0.075329567, -0.071428571, -0.063989962, -0.051378446, -0.03081761, 0];
Vlo_8 = Ulo_8.*Cyl;
Ulo_8 = Ulo_8.*Cxl;

Ulo_10 = [-0.475015813, 0.661587302, 0.371900826, 0.225602028, 0.116899619, 0.049618321, 0.004444444, -0.039974619, -0.060991105, -0.079415502, -0.091486658, -0.105463787, -0.116190476, -0.125634518, -0.137142857, -0.146295123, -0.142857143, -0.146853147, 0];
Vlo_10 = Ulo_10.*Cyl;
Ulo_10 = Ulo_10.*Cxl;

Ulo_12 = [-0.000632511, 0.625631313, 0.335863378, 0.193303853, 0.081115336, 0.003792668, -0.044303797, -0.092346616, -0.124368687, -0.142586751, -0.163084703, -0.185958254, -0.207726409, -0.224257738, -0.248894504, -0.270714738, -0.292035398, -0.310126582, 0];
Vlo_12 = Ulo_12.*Cyl;
Ulo_12 = Ulo_12.*Cxl;

Ulo_14 = [0.021410579, 0.644612476, 0.353943218, 0.207939509, 0.090334807, 0.00754717, -0.044164038, -0.094065657, -0.111111111, -0.154040404, -0.18147448, -0.208832808, -0.230429293, -0.255359395, -0.282471627, -0.309343434, -0.338587642, -0.355387524, 0];
Vlo_14 = Ulo_14.*Cyl;
Ulo_14 = Ulo_14.*Cxl;

% chordwise and chordnormal graph calcs
Uup_0_cw = Cwx.*Uup_0c;
Vup_0_cw = Cwy.*Vup_0c;
Uup_0_cn = Cnx.*Uup_0c;
Vup_0_cn = Cny.*Vup_0c;
Ulo_0_cw = Cwx.*Ulo_0;
Vlo_0_cw = Cwy.*Vlo_0;
Ulo_0_cn = Cnx.*Ulo_0;
Vlo_0_cn = Cny.*Vlo_0;

Uup_8_cw = Cwx.*Uup_8;
Vup_8_cw = Cwy.*Vup_8;
Uup_8_cn = Cnx.*Uup_8;
Vup_8_cn = Cny.*Vup_8;
Ulo_8_cw = Cwx.*Ulo_8;
Vlo_8_cw = Cwy.*Vlo_8;
Ulo_8_cn = Cnx.*Ulo_8;
Vlo_8_cn = Cny.*Vlo_8;

Uup_10_cw = Cwx.*Uup_10;
Vup_10_cw = Cwy.*Vup_10;
Uup_10_cn = Cnx.*Uup_10;
Vup_10_cn = Cny.*Vup_10;
Ulo_10_cw = Cwx.*Ulo_10;
Vlo_10_cw = Cwy.*Vlo_10;
Ulo_10_cn = Cnx.*Ulo_10;
Vlo_10_cn = Cny.*Vlo_10;

Uup_12_cw = Cwx.*Uup_12;
Vup_12_cw = Cwy.*Vup_12;
Uup_12_cn = Cnx.*Uup_12;
Vup_12_cn = Cny.*Vup_12;
Ulo_12_cw = Cwx.*Ulo_12;
Vlo_12_cw = Cwy.*Vlo_12;
Ulo_12_cn = Cnx.*Ulo_12;
Vlo_12_cn = Cny.*Vlo_12;

% plots
% 10deg aoa cp vectors
figure(1);
hold on
axis equal
quiver(Xup, Yup, Uup_10, Vup_10, "color", "black");
quiver(Xlo, Ylo, Ulo_10, Vlo_10, "color", "black");
plot(Xup, Yup, "color", "black");
plot(Xlo, Ylo, "color", "black");
title("10 deg AoA Cp")
hold off
% 0deg aoa cp cw vectors
figure(2);
hold on
axis equal
quiver(Xup, Yup, Uup_0_cw, Vup_0_cw, "color", "black");
quiver(Xlo, Ylo, Ulo_0_cw, Vlo_0_cw, "color", "black");
plot(Xup, Yup, "color", "black");
plot(Xlo, Ylo, "color", "black");
title("0 deg Aoa Cp Chordwise")
hold off
% 0 deg aoa cp cn vectors
figure(3);
hold on
axis equal
quiver(Xup, Yup, Uup_0_cn, Vup_0_cn, "color", "black");
quiver(Xlo, Ylo, Ulo_0_cn, Vlo_0_cn, "color", "black");
plot(Xup, Yup, "color", "black");
plot(Xlo, Ylo, "color", "black");
title("0 deg Aoa Cp Chord Normal")
hold off
% 8 deg aoa cp cw vectors
figure(4);
hold on
axis equal
quiver(Xup, Yup, Uup_8_cw, Vup_8_cw, "color", "black");
quiver(Xlo, Ylo, Ulo_8_cw, Vlo_8_cw, "color", "black");
plot(Xup, Yup, "color", "black");
plot(Xlo, Ylo, "color", "black");
title("8 deg AoA Cp Chordwise");
hold off
% 8 deg aoa cp cn vectors
figure(5);
hold on
axis equal
quiver(Xup, Yup, Uup_8_cn, Vup_8_cn, "color", "black");
quiver(Xlo, Ylo, Ulo_8_cn, Vlo_8_cn, "color", "black");
plot(Xup, Yup, "color", "black");
plot(Xlo, Ylo, "color", "black");
title("8 deg AoA Cp Chord Normal");
hold off
% 10 deg aoa cp cw vectors
figure(6);
hold on
axis equal
quiver(Xup, Yup, Uup_10_cw, Vup_10_cw, "color", "black");
quiver(Xlo, Ylo, Ulo_10_cw, Vlo_10_cw, "color", "black");
plot(Xup, Yup, "color", "black");
plot(Xlo, Ylo, "color", "black");
title("10 deg AoA Cp Chordwise")
hold off
% 10 deg aoa cp cn vectors
figure(7);
hold on
axis equal
quiver(Xup, Yup, Uup_10_cn, Vup_10_cn, "color", "black");
quiver(Xlo, Ylo, Ulo_10_cn, Vlo_10_cn, "color", "black");
plot(Xup, Yup, "color", "black");
plot(Xlo, Ylo, "color", "black");
title("10 deg AoA Cp Chord Normal")
hold off
% 12 deg aoa cp cw vectors
figure(8);
hold on
axis equal
quiver(Xup, Yup, Uup_12_cw, Vup_12_cw, "color", "black");
quiver(Xlo, Ylo, Ulo_12_cw, Vlo_12_cw, "color", "black");
plot(Xup, Yup, "color", "black");
plot(Xlo, Ylo, "color", "black");
title("12 deg AoA Cp Chordwise")
hold off
% 12 deg aoa cp cn vectors
figure(9);
hold on
axis equal
quiver(Xup, Yup, Uup_12_cn, Vup_12_cn, "color", "black");
quiver(Xlo, Ylo, Ulo_12_cn, Vlo_12_cn, "color", "black");
plot(Xup, Yup, "color", "black");
plot(Xlo, Ylo, "color", "black");
title("12 deg AoA Cp Chord Normal")
hold off