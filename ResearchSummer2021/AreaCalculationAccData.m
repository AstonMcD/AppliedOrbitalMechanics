t = [-1.509590983390808, -1.409936010837555, -1.3102799654006958, -1.2106240391731262, -1.1109690070152283, -1.0113129615783691, -0.9116570353507996, -0.8120030164718628, -0.7123469710350037, -0.6126909852027893, -0.5130349397659302, -0.413379967212677, -0.3137240409851074, -0.21406793594360352, -0.11441296339035034, -0.014757037162780762, 0.08489900827407837, 0.18455404043197632, 0.2842099666595459, 0.38386601209640503, 0.4835200309753418];
linx = [13.273195785134964, 8.027501679655307, 3.1792085995724753, 11.206710234165628, 5.961016109522698, 9.855546600672568, 6.199456752613948, 12.6373540779612, 267.0534602556749, 675.1840271445798, 304.1707020857225, -83.61319257560625, -43.47567938973686, 21.539137913240125, 8.504382963865913, 0.0, 2.5433668811026795, 14.465398984118115, 4.132971175982902, 3.1792085995724753, 4.371411819799212];
plot(t,linx)
% title("Change in Acceleration with time for THR_0")
% xlabel("Time")
% ylabel("Linear x-Acceleration")
Area = trapz(t,linx)