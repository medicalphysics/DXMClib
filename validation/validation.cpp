

#include "dxmc/tube.h"
#include "dxmc/source.h"
#include "dxmc/world.h"
#include "dxmc/transport.h"

#include <numeric>
#include <iostream>
#include <cassert>

constexpr double ERRF = 1e-4;


// energy weighs pair for spectre
/*RQR-8
W/Al
100 kVp
11 deg anode angle
0% Ripple
Al filter thickness: 2.708 mm
Mean Energy: 50.6 keV
HVL: 3.950 mm Al
QVL: 9.840 mm Al
*/
const std::vector<double> TG195_100KV({ 16.25,1.423E-04,16.75,2.157E-04,17.25,3.102E-04,17.75,4.324E-04,18.25,5.840E-04,18.75,7.644E-04,19.25,9.784E-04,19.75,1.222E-03,20.25,1.491E-03,20.75,1.803E-03,21.25,2.129E-03,21.75,2.490E-03,22.25,2.863E-03,22.75,3.263E-03,23.25,3.658E-03,23.75,4.093E-03,24.25,4.504E-03,24.75,4.912E-03,25.25,5.347E-03,25.75,5.769E-03,26.25,6.168E-03,26.75,6.582E-03,27.25,6.965E-03,27.75,7.360E-03,28.25,7.710E-03,28.75,8.067E-03,29.25,8.368E-03,29.75,8.671E-03,30.25,8.975E-03,30.75,9.213E-03,31.25,9.476E-03,31.75,9.694E-03,32.25,9.903E-03,32.75,1.009E-02,33.25,1.025E-02,33.75,1.040E-02,34.25,1.053E-02,34.75,1.063E-02,35.25,1.073E-02,35.75,1.081E-02,36.25,1.087E-02,36.75,1.092E-02,37.25,1.096E-02,37.75,1.099E-02,38.25,1.100E-02,38.75,1.100E-02,39.25,1.099E-02,39.75,1.098E-02,40.25,1.095E-02,40.75,1.091E-02,41.25,1.086E-02,41.75,1.081E-02,42.25,1.076E-02,42.75,1.069E-02,43.25,1.063E-02,43.75,1.055E-02,44.25,1.048E-02,44.75,1.039E-02,45.25,1.031E-02,45.75,1.022E-02,46.25,1.012E-02,46.75,1.003E-02,47.25,9.933E-03,47.75,9.828E-03,48.25,9.732E-03,48.75,9.628E-03,49.25,9.516E-03,49.75,9.412E-03,50.25,9.302E-03,50.75,9.193E-03,51.25,9.084E-03,51.75,8.970E-03,52.25,8.862E-03,52.75,8.749E-03,53.25,8.637E-03,53.75,8.526E-03,54.25,8.409E-03,54.75,8.300E-03,55.25,8.185E-03,55.75,8.072E-03,56.25,7.959E-03,56.75,7.847E-03,57.25,7.737E-03,57.75,2.568E-02,58.25,7.513E-03,58.75,7.405E-03,59.25,3.920E-02,59.75,7.181E-03,60.25,7.071E-03,60.75,6.962E-03,61.25,6.854E-03,61.75,6.746E-03,62.25,6.640E-03,62.75,6.530E-03,63.25,6.425E-03,63.75,6.321E-03,64.25,6.214E-03,64.75,6.107E-03,65.25,6.006E-03,65.75,5.901E-03,66.25,5.797E-03,66.75,1.673E-02,67.25,5.592E-03,67.75,5.491E-03,68.25,5.390E-03,68.75,8.223E-03,69.25,5.055E-03,69.75,4.296E-03,70.25,4.236E-03,70.75,4.171E-03,71.25,4.110E-03,71.75,4.048E-03,72.25,3.982E-03,72.75,3.919E-03,73.25,3.852E-03,73.75,3.787E-03,74.25,3.719E-03,74.75,3.654E-03,75.25,3.585E-03,75.75,3.516E-03,76.25,3.449E-03,76.75,3.379E-03,77.25,3.308E-03,77.75,3.240E-03,78.25,3.169E-03,78.75,3.098E-03,79.25,3.026E-03,79.75,2.954E-03,80.25,2.882E-03,80.75,2.809E-03,81.25,2.736E-03,81.75,2.665E-03,82.25,2.592E-03,82.75,2.519E-03,83.25,2.445E-03,83.75,2.370E-03,84.25,2.296E-03,84.75,2.222E-03,85.25,2.148E-03,85.75,2.073E-03,86.25,1.999E-03,86.75,1.925E-03,87.25,1.850E-03,87.75,1.776E-03,88.25,1.700E-03,88.75,1.625E-03,89.25,1.550E-03,89.75,1.476E-03,90.25,1.400E-03,90.75,1.326E-03,91.25,1.251E-03,91.75,1.177E-03,92.25,1.101E-03,92.75,1.027E-03,93.25,9.529E-04,93.75,8.781E-04,94.25,8.041E-04,94.75,7.302E-04,95.25,6.559E-04,95.75,5.823E-04,96.25,5.089E-04,96.75,4.353E-04,97.25,3.623E-04,97.75,2.892E-04,98.25,2.166E-04,98.75,1.441E-04,99.25,7.193E-05,99.75,5.990E-06 });

/*RQR-9
W/Al
120 kVp
11 deg anode angle
0% Ripple
Al filter thickness: 2.861 mm
Mean Energy: 56.4 keV
HVL: 5.010 mm Al
*/
const std::vector<double> TG195_120KV({ 16.75,1.107E-04,17.25,1.625E-04,17.75,2.308E-04,18.25,3.172E-04,18.75,4.220E-04,19.25,5.486E-04,19.75,6.956E-04,20.25,8.610E-04,20.75,1.056E-03,21.25,1.264E-03,21.75,1.499E-03,22.25,1.748E-03,22.75,2.019E-03,23.25,2.293E-03,23.75,2.601E-03,24.25,2.900E-03,24.75,3.203E-03,25.25,3.531E-03,25.75,3.858E-03,26.25,4.176E-03,26.75,4.511E-03,27.25,4.830E-03,27.75,5.163E-03,28.25,5.469E-03,28.75,5.786E-03,29.25,6.065E-03,29.75,6.349E-03,30.25,6.638E-03,30.75,6.879E-03,31.25,7.143E-03,31.75,7.372E-03,32.25,7.597E-03,32.75,7.804E-03,33.25,7.994E-03,33.75,8.171E-03,34.25,8.339E-03,34.75,8.483E-03,35.25,8.622E-03,35.75,8.745E-03,36.25,8.849E-03,36.75,8.949E-03,37.25,9.031E-03,37.75,9.109E-03,38.25,9.170E-03,38.75,9.219E-03,39.25,9.264E-03,39.75,9.297E-03,40.25,9.319E-03,40.75,9.332E-03,41.25,9.333E-03,41.75,9.332E-03,42.25,9.327E-03,42.75,9.307E-03,43.25,9.292E-03,43.75,9.259E-03,44.25,9.229E-03,44.75,9.187E-03,45.25,9.149E-03,45.75,9.101E-03,46.25,9.044E-03,46.75,8.996E-03,47.25,8.937E-03,47.75,8.871E-03,48.25,8.813E-03,48.75,8.747E-03,49.25,8.672E-03,49.75,8.605E-03,50.25,8.530E-03,50.75,8.456E-03,51.25,8.381E-03,51.75,8.300E-03,52.25,8.226E-03,52.75,8.145E-03,53.25,8.065E-03,53.75,7.985E-03,54.25,7.899E-03,54.75,7.820E-03,55.25,7.736E-03,55.75,7.652E-03,56.25,7.568E-03,56.75,7.486E-03,57.25,7.403E-03,57.75,3.335E-02,58.25,7.236E-03,58.75,7.155E-03,59.25,5.339E-02,59.75,6.986E-03,60.25,6.903E-03,60.75,6.821E-03,61.25,6.739E-03,61.75,6.658E-03,62.25,6.578E-03,62.75,6.494E-03,63.25,6.415E-03,63.75,6.338E-03,64.25,6.256E-03,64.75,6.175E-03,65.25,6.100E-03,65.75,6.021E-03,66.25,5.942E-03,66.75,2.242E-02,67.25,5.788E-03,67.75,5.712E-03,68.25,5.637E-03,68.75,9.988E-03,69.25,5.257E-03,69.75,4.045E-03,70.25,4.019E-03,70.75,3.988E-03,71.25,3.960E-03,71.75,3.932E-03,72.25,3.900E-03,72.75,3.871E-03,73.25,3.838E-03,73.75,3.808E-03,74.25,3.774E-03,74.75,3.743E-03,75.25,3.709E-03,75.75,3.674E-03,76.25,3.641E-03,76.75,3.606E-03,77.25,3.570E-03,77.75,3.537E-03,78.25,3.500E-03,78.75,3.463E-03,79.25,3.426E-03,79.75,3.389E-03,80.25,3.351E-03,80.75,3.313E-03,81.25,3.274E-03,81.75,3.238E-03,82.25,3.200E-03,82.75,3.160E-03,83.25,3.121E-03,83.75,3.079E-03,84.25,3.039E-03,84.75,3.000E-03,85.25,2.959E-03,85.75,2.919E-03,86.25,2.878E-03,86.75,2.838E-03,87.25,2.797E-03,87.75,2.756E-03,88.25,2.712E-03,88.75,2.671E-03,89.25,2.629E-03,89.75,2.588E-03,90.25,2.544E-03,90.75,2.502E-03,91.25,2.460E-03,91.75,2.418E-03,92.25,2.374E-03,92.75,2.331E-03,93.25,2.289E-03,93.75,2.244E-03,94.25,2.202E-03,94.75,2.159E-03,95.25,2.115E-03,95.75,2.072E-03,96.25,2.029E-03,96.75,1.984E-03,97.25,1.941E-03,97.75,1.896E-03,98.25,1.853E-03,98.75,1.809E-03,99.25,1.765E-03,99.75,1.722E-03,100.25,1.677E-03,100.75,1.634E-03,101.25,1.589E-03,101.75,1.546E-03,102.25,1.501E-03,102.75,1.458E-03,103.25,1.414E-03,103.75,1.370E-03,104.25,1.326E-03,104.75,1.282E-03,105.25,1.238E-03,105.75,1.195E-03,106.25,1.151E-03,106.75,1.107E-03,107.25,1.063E-03,107.75,1.019E-03,108.25,9.761E-04,108.75,9.323E-04,109.25,8.893E-04,109.75,8.456E-04,110.25,8.027E-04,110.75,7.592E-04,111.25,7.158E-04,111.75,6.731E-04,112.25,6.300E-04,112.75,5.874E-04,113.25,5.445E-04,113.75,5.017E-04,114.25,4.594E-04,114.75,4.168E-04,115.25,3.747E-04,115.75,3.324E-04,116.25,2.903E-04,116.75,2.485E-04,117.25,2.067E-04,117.75,1.650E-04,118.25,1.236E-04,118.75,8.222E-05,119.25,4.102E-05,119.75,3.417E-06 });



bool inline isEqual(double a, double b)
{
	return std::abs(a - b) < ERRF;
}


World generateTG195Case2World()
{

	Material air("C0.02N78.44O21.08Ar0.47");
	air.setStandardDensity(0.001205);
	Material soft("H62.95C12.91N1.17O22.78Na0.03P0.04S0.06Cl0.03K0.03");
	soft.setStandardDensity(1.03);

	std::array<double, 3> spacing = { 1,1,200 }; // mm
	std::array<std::size_t, 3> dim = { 390,390,9 };
	const auto size = dim[0] * dim[1] * dim[2];
	auto dens = std::make_shared<std::vector<double>>(size, air.standardDensity());
	auto idx = std::make_shared<std::vector<unsigned char>>(size, 0);
	// fill first slice with soft tissue
	for (std::size_t i = 0; i < dim[0] * dim[1]; ++i)
	{
		dens->data()[i] = soft.standardDensity();
		idx->data()[i] = static_cast<unsigned char>(1);
	}


	World w;
	w.setDimensions(dim);
	w.setSpacing(spacing);
	w.setDensityArray(dens);
	w.addMaterialToMap(air);
	w.addMaterialToMap(soft);
	w.setMaterialIndexArray(idx);
	return w;
}


bool TG195Case2AbsorbedEnergyMono()
{
	auto w = generateTG195Case2World();
	IsotropicSource src;
	std::vector<double> s({ 1.0 }), e({ 56.4 });
	w.setAttenuationLutMaxEnergy(60.0);
	//w.setAttenuationLutResolution(0.4);
	w.validate();
	assert(w.isValid());


	src.setSpecter(s, e);
	auto extent = w.matrixExtent();
	src.setPosition(0, 0, 1750 + extent[4]);
	const double halfAng = std::atan(390 * 0.5 / 1800.0);
	src.setCollimationAngles(halfAng * 2.0, halfAng * 2.0);
	std::array<double, 6> cosines = { -1,0,0,0,1,0 };
	src.setDirectionCosines(cosines);
	src.setHistoriesPerExposure(1e6);
	src.setTotalExposures(16);
	src.validate();
	auto dose = transport::run(w, &src, nullptr, false);

	const double total_kev = std::accumulate(dose.cbegin(), dose.cend(), 0.0);
	const double total_hist = static_cast<double>(src.totalExposures() * src.historiesPerExposure());

	const double ev_per_hist = total_kev * 1000.0 / total_hist;
	const double sim_ev_per_hist = 33171.4;
	std::cout << "TG195 Case 2\n";
	std::cout << "Monochromatic source of 56.4 keV:\n";
	std::cout << "eV per history: " << ev_per_hist << "   real eV per history: " << sim_ev_per_hist << "\n";
	std::cout << "deviation: " << ev_per_hist- sim_ev_per_hist << " (%): " << (ev_per_hist - sim_ev_per_hist)/sim_ev_per_hist*100 << "\n";

	return true;
}

bool TG195Case2AbsorbedEnergy120()
{
	auto w = generateTG195Case2World();
	IsotropicSource src;

	std::vector<double> s, e;
	for (std::size_t i = 0; i < TG195_120KV.size(); i = i + 2)
	{
		e.push_back(TG195_120KV[i]);
		s.push_back(TG195_120KV[i + 1]);
	}
	const auto max_energy = *(std::max_element(e.cbegin(), e.cend()));
	w.setAttenuationLutMaxEnergy(max_energy);
	//w.setAttenuationLutResolution(0.4);
	w.validate();
	assert(w.isValid());
	src.setSpecter(s, e);
	auto extent = w.matrixExtent();
	src.setPosition(0, 0, 1750 + extent[4]);
	const double halfAng = std::atan(390 * 0.5 / 1800.0);
	src.setCollimationAngles(halfAng * 2.0, halfAng * 2.0);
	std::array<double, 6> cosines = { -1,0,0,0,1,0 };
	src.setDirectionCosines(cosines);
	src.setHistoriesPerExposure(1e6);
	src.setTotalExposures(16);
	src.validate();
	auto dose = transport::run(w, &src, nullptr, false);

	const double total_kev = std::accumulate(dose.cbegin(), dose.cend(), 0.0);
	const double total_hist = static_cast<double>(src.totalExposures() * src.historiesPerExposure());

	const double ev_per_hist = total_kev * 1000.0 / total_hist;
	const double sim_ev_per_hist = 33125.98;
	std::cout << "TG195 Case 2\n";
	std::cout << "Specter source of 120 kV:\n";
	std::cout << "eV per history: " << ev_per_hist << "   real eV per history: " << sim_ev_per_hist << "\n";
	std::cout << "deviation: " << ev_per_hist - sim_ev_per_hist << " (%): " << (ev_per_hist - sim_ev_per_hist) / sim_ev_per_hist * 100 << "\n";

	return true;
}




int main(int argc, char* argv[])
{

	auto success = TG195Case2AbsorbedEnergy120();
	success = success && TG195Case2AbsorbedEnergyMono();
	if (success)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}