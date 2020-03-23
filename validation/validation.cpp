

#include "dxmc/tube.h"
#include "dxmc/source.h"
#include "dxmc/world.h"
#include "dxmc/transport.h"

#include <numeric>
#include <execution>
#include <iostream>
#include <cassert>
#include <fstream>

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

bool betw(std::size_t v, std::size_t min, std::size_t max)
{
	if (v >= min && v < max)
		return true;
	return false;
}
bool betw(double v, double min, double max)
{
	if (v >= min && v < max)
		return true;
	return false;
}
World generateTG195Case2World()
{

	Material air("C0.02N78.44O21.08Ar0.47");
	air.setStandardDensity(0.001205);
	Material soft("H62.95C12.91N1.17O22.78Na0.03P0.04S0.06Cl0.03K0.03");
	soft.setStandardDensity(1.03);

	std::array<double, 3> spacing = { 5,5,5 }; // mm
	std::array<std::size_t, 3> dim = { 78,78,400};
	const auto size = dim[0] * dim[1] * dim[2];
	auto dens = std::make_shared<std::vector<double>>(size, air.standardDensity());
	auto idx = std::make_shared<std::vector<unsigned char>>(size, 0);
	// fill first slice with soft tissue
	for (std::size_t z = 0; z < dim[2]; ++z)
		for (std::size_t y = 0; y < dim[1]; ++y)
			for (std::size_t x = 0; x < dim[0]; ++x) {
				if (z < 40)
				{
					const auto i = x + y * dim[0] + z * dim[0] * dim[1];
					dens->data()[i] = soft.standardDensity();
					idx->data()[i] = static_cast<unsigned char>(1);

					//center boxes
					if (betw(x, 36, 42) && betw(y, 36, 42) && betw(z, 5, 11))
						idx->data()[i] = static_cast<unsigned char>(9+1);
					if (betw(x, 36, 42) && betw(y, 36, 42) && betw(z, 11, 17))
						idx->data()[i] = static_cast<unsigned char>(8+1);
					if (betw(x, 36, 42) && betw(y, 36, 42) && betw(z, 17, 23))
						idx->data()[i] = static_cast<unsigned char>(3+1);
					if (betw(x, 36, 42) && betw(y, 36, 42) && betw(z, 23, 29))
						idx->data()[i] = static_cast<unsigned char>(7+1);
					if (betw(x, 36, 42) && betw(y, 36, 42) && betw(z, 29, 35))
						idx->data()[i] = static_cast<unsigned char>(6+1);
					
					//periphery y
					if (betw(x, 36, 42) && betw(y, 6, 12) && betw(z, 17, 23))
						idx->data()[i] = static_cast<unsigned char>(1+1);
					if (betw(x, 36, 42) && betw(y, 66, 72) && betw(z, 17, 23))
						idx->data()[i] = static_cast<unsigned char>(5 + 1);
					if (betw(x, 6, 12) && betw(y, 36, 42) && betw(z, 17, 23))
						idx->data()[i] = static_cast<unsigned char>(2 + 1);
					if (betw(x, 66, 72) && betw(y, 36, 42) && betw(z, 17, 23))
						idx->data()[i] = static_cast<unsigned char>(4 + 1);
				}

			}
	/*std::array<double, 3> cm = { 0,0,0 };
	double m = 0;
	for (std::size_t z = 0; z < dim[2]; ++z)
		for (std::size_t y = 0; y < dim[1]; ++y)
			for (std::size_t x = 0; x < dim[0]; ++x) {
				const auto i = x + y * dim[0] + z * dim[0] * dim[1];
				if (idx->data()[i] == 5 + 1)
				{
					m += 1.0;
					cm[0] += spacing[0] * x;
					cm[1] += spacing[1] * y;
					cm[2] += spacing[2] * z;
				}
			}
	std::cout << "mass center: " << cm[0] / m << ", " << cm[1] / m << ", "<<cm[2] / m << "\n";
	*/
	World w;
	w.setDimensions(dim);
	w.setSpacing(spacing);
	w.setDensityArray(dens);
	w.addMaterialToMap(air);
	for (std::size_t i = 0; i < 10; ++i)
		w.addMaterialToMap(soft);
	w.setMaterialIndexArray(idx);
	return w;
}


bool TG195Case2AbsorbedEnergyMono()
{
	std::cout << "TG195 Case 2\n";
	std::cout << "Monochromatic source of 56.4 keV:\n";
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
	src.setTotalExposures(40);
	src.validate();
	auto dose = transport::run(w, &src, nullptr, false);

	const double total_hist = static_cast<double>(src.totalExposures() * src.historiesPerExposure());
	double total_ev = 1000*std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [](double d, unsigned char i)->double {return i > 0 ? d : 0.0; });
	total_ev /= total_hist;
	std::array<double, 9> subvol_ev;
	for (std::size_t i = 0; i < 9; ++i) {
		subvol_ev[i] = 1000*std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m)->double {return m == i + 2 ? d : 0.0; });
		subvol_ev[i] /= total_hist;
	}

	const double sim_ev = 33171.4;
	std::array<double, 9> sim_subvol = { 27.01,27.00,36.67,27.01,27.01,72.86,53.35,23.83,14.60 };

	std::cout << "Total body dxmc: " << total_ev << ", TG195: " << sim_ev << ", difference [%]: " << (total_ev - sim_ev) / sim_ev * 100 << "\n";
	for (std::size_t i =0;i < subvol_ev.size();++i)
		std::cout << "VOI "<< i+1 << " dxmc: " << subvol_ev[i] << ", TG195: " << sim_subvol[i] << ", difference [%]: " << (subvol_ev[i] - sim_subvol[i]) / sim_subvol[i] * 100 << "\n";
	std::cout << "\n";
	return true;
}

bool TG195Case2AbsorbedEnergy120()
{
	std::cout << "TG195 Case 2\n";
	std::cout << "Specter source of 120 kV:\n";

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
	src.setTotalExposures(40);
	src.validate();
	auto dose = transport::run(w, &src, nullptr, false);

	const double total_hist = static_cast<double>(src.totalExposures() * src.historiesPerExposure());
	double total_ev = 1000 * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [](double d, unsigned char i)->double {return i > 0 ? d : 0.0; });
	total_ev /= total_hist;
	std::array<double, 9> subvol_ev;
	for (std::size_t i = 0; i < 9; ++i) {
		subvol_ev[i] = 1000 * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m)->double {return m == i + 2 ? d : 0.0; });
		subvol_ev[i] /= total_hist;
	}

	const double sim_ev = 33125.98;
	std::array<double, 9> sim_subvol = { 24.97,24.95,33.52,24.96,24.97,72.70,49.99,21.73,13.48};

	std::cout << "Total body dxmc: " << total_ev << ", TG195: " << sim_ev << ", difference [%]: " << (total_ev - sim_ev) / sim_ev * 100 << "\n";
	for (std::size_t i = 0; i < subvol_ev.size(); ++i)
		std::cout << "VOI " << i + 1 << " dxmc: " << subvol_ev[i] << ", TG195: " << sim_subvol[i] << ", difference [%]: " << (subvol_ev[i] - sim_subvol[i]) / sim_subvol[i] * 100 << "\n";
	std::cout << "\n";
	return true;
}

std::vector<std::size_t> circleIndices(const double center_x, const double center_y, const std::array<std::size_t, 3>& dim, const std::array<double, 3>& spacing, const double radii=1.0)
{
	std::vector<std::size_t> ind;
	const double radii2 = radii * radii;
	for (std::size_t yi = 0; yi < dim[1]; ++yi)
		for (std::size_t xi = 0; xi < dim[0]; ++xi)
		{
			const double x = xi * spacing[0] - center_x;
			const double y = yi * spacing[1] - center_y;
			if ((x*x+y*y) <= radii2)
			{
				const auto i = xi  + yi * dim[0];
				ind.push_back(i);
			}
		}
	return ind;
}

World generateTG195Case4World1() {

	const std::array<std::size_t, 3> dim = {320,320,600};
	const std::array<double, 3> spacing = { 1,1,5 };
	const auto size = std::accumulate(dim.cbegin(), dim.cend(), 1.0, std::multiplies<>());

	Material air("C0.02N78.44O21.08Ar0.47");
	air.setStandardDensity(0.001205);
	Material pmma("H53.28C33.37O13.35");
	pmma.setStandardDensity(1.19);

	auto mat = std::make_shared<std::vector<unsigned char>>(size, static_cast<unsigned char>(0));
	auto dens = std::make_shared<std::vector<double>>(size, air.standardDensity());
	//generate cylindar
	auto circ_ind = circleIndices(160.0, 160.0, dim, spacing, 160.0);
	for (std::size_t z =0;z < dim[2];++z)
		for (const auto ind : circ_ind)
		{
			const auto i = z * dim[0] * dim[1] + ind;
			(mat->data())[i] = static_cast<unsigned char>(1);
			(dens->data())[i] = pmma.standardDensity();
			if (betw((z + 0.5) * spacing[2] - spacing[2] * dim[2] * 0.5, -5, 5))
				mat->data()[i] = static_cast<unsigned char>(2);
			if (betw((z + 0.5) * spacing[2] - spacing[2] * dim[2] * 0.5, -15, -5))
				mat->data()[i] = static_cast<unsigned char>(3);
			if (betw((z + 0.5) * spacing[2] - spacing[2] * dim[2] * 0.5, -25, -15))
				mat->data()[i] = static_cast<unsigned char>(4);
			if (betw((z + 0.5) * spacing[2] - spacing[2] * dim[2] * 0.5, -35, -25))
				mat->data()[i] = static_cast<unsigned char>(5);
		}


	World w;
	w.setSpacing(spacing);
	w.setDimensions(dim);
	w.setDensityArray(dens);
	w.setMaterialIndexArray(mat);
	w.addMaterialToMap(air);
	w.addMaterialToMap(pmma);
	for (int i = 0; i < 4; ++i)
		w.addMaterialToMap(pmma);
	w.validate();
	return w;
}

World generateTG195Case4World2() {

	const std::array<std::size_t, 3> dim = { 320,320,600 };
	const std::array<double, 3> spacing = { 1,1,5 };
	const auto size = std::accumulate(dim.cbegin(), dim.cend(), 1.0, std::multiplies<>());

	Material air("C0.02N78.44O21.08Ar0.47");
	air.setStandardDensity(0.001205);
	Material pmma("H53.28C33.37O13.35");
	pmma.setStandardDensity(1.19);

	auto mat = std::make_shared<std::vector<unsigned char>>(size, static_cast<unsigned char>(0));
	auto dens = std::make_shared<std::vector<double>>(size, air.standardDensity());
	//generate cylindar
	auto circ_ind = circleIndices(320. * .5, 320. * .5, dim, spacing, 320. * 0.5);
	auto pmma2_ind = circleIndices(320. * .5 - 150, 320. * .5, dim, spacing, 10.0);
	auto pmma1_ind = circleIndices(320. * .5, 320. * .5, dim, spacing, 10.0);
	for (std::size_t z = 0; z < dim[2]; ++z)
		for (const auto i : circ_ind)
		{
			const auto ind = z * dim[0] * dim[1] + i;
			mat->data()[ind] = static_cast<unsigned char>(1);
			dens->data()[ind] = pmma.standardDensity();
		}
	for (const auto i : pmma1_ind)
	{
		const auto ind1 = 299*  dim[0] * dim[1] + i;
		const auto ind2 = 300 * dim[0] * dim[1] + i;
		mat->data()[ind1] = static_cast<unsigned char>(2);
		mat->data()[ind2] = static_cast<unsigned char>(2);
	}

	for (const auto i : pmma2_ind)
	{
		const auto ind1 = 299 * dim[0] * dim[1] + i;
		const auto ind2 = 300 * dim[0] * dim[1] + i;
		mat->data()[ind1] = static_cast<unsigned char>(3);
		mat->data()[ind2] = static_cast<unsigned char>(3);
	}

	World w;
	w.setSpacing(spacing);
	w.setDimensions(dim);
	w.setDensityArray(dens);
	w.setMaterialIndexArray(mat);
	w.addMaterialToMap(air);
	w.addMaterialToMap(pmma);
	for (int i = 0; i < 2; ++i)
		w.addMaterialToMap(pmma);
	w.validate();
	return w;
}


bool TG195Case4AbsorbedEnergy() {
	IsotropicSource src;
	src.setPosition(-600.0, 0., 0.);
	std::array<double, 6> cos = { 0,1,0 , 0,0,1 };
	src.setDirectionCosines(cos);
	
	double direction[3];
	vectormath::cross(cos.data(), direction);

	std::cout << "TG195 Case 4.1:\nMonoenergetic specter of 56.4 kev\n";
	std::cout << "Collimation: 10 mm:\n";
	std::vector<double> s, e;
	s.push_back(1.0);
	e.push_back(56.4);
	src.setSpecter(s, e);
	src.setHistoriesPerExposure(1e6);
	src.setTotalExposures(16);
	src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(5. / 600.) * 2.);
	auto w = generateTG195Case4World1();

	auto dose = transport::run(w, &src, nullptr, false);

	const double total_hist = static_cast<double>(src.totalExposures() * src.historiesPerExposure());
	std::array<double, 4> voi_ev;
	for (std::size_t i = 0; i < 4; ++i) {
		voi_ev[i] = 1000.0 / total_hist * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m)->double {return m == i + 2 ? d : 0.0; });
	}
	std::array<double, 4> sim_ev = { 11592.27,2576.72,1766.85,1330.53 };
	for (std::size_t i = 0; i < voi_ev.size(); ++i)
		std::cout << "VOI " << i + 1 << " dxmc: " << voi_ev[i] << ", TG195: " << sim_ev[i] << ", difference [%]: " << (voi_ev[i] - sim_ev[i]) / sim_ev[i] * 100 << "\n";
	std::cout << std::endl;
	
	
	/*auto myfile = std::fstream("dens.bin", std::ios::out | std::ios::binary | std::ios::trunc);
	myfile.write((const char*)(w.densityBuffer()), sizeof(double) * dose.size());
	myfile.close();
	myfile = std::fstream("mat.bin", std::ios::out | std::ios::binary | std::ios::trunc);
	myfile.write((const char*)(w.materialIndexBuffer()), sizeof(unsigned char) * dose.size());
	myfile.close();
	myfile = std::fstream("dose.bin", std::ios::out | std::ios::binary | std::ios::trunc);
	myfile.write((const char*)(dose.data()), sizeof(double) * dose.size());
	myfile.close();*/


	std::cout << "Collimation: 80 mm:\n";
	src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(40. / 600.) * 2.);
	
	dose = transport::run(w, &src, nullptr, false);

	for (std::size_t i = 0; i < 4; ++i) {
		voi_ev[i] = 1000 / total_hist * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m)->double {return m == i + 2 ? d : 0.0; });
	}

	sim_ev = { 3380.39,3332.64,3176.44,2559.58 };

	for (std::size_t i = 0; i < voi_ev.size(); ++i)
		std::cout << "VOI " << i + 1 << " dxmc: " << voi_ev[i] << ", TG195: " << sim_ev[i] << ", difference [%]: " << (voi_ev[i] - sim_ev[i]) / sim_ev[i] * 100 << "\n";



	std::cout << "TG195 Case 4.1:\nSpecter of 120 kV W/Al\n";
	std::cout << "Collimation: 10 mm:\n";
	   	  
	for (std::size_t i = 0; i < TG195_120KV.size(); i = i + 2)
	{
		e.push_back(TG195_120KV[i]);
		s.push_back(TG195_120KV[i + 1]);
	}
	src.setSpecter(s, e);
	src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(5. / 600.) * 2.);
	dose = transport::run(w, &src, nullptr, false);
	for (std::size_t i = 0; i < 4; ++i) {
		voi_ev[i] = 1000 / total_hist * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m)->double {return m == i + 2 ? d : 0.0; });
	}
	sim_ev = { 13137.02,2585.47,1706.86,1250.61	};

	for (std::size_t i = 0; i < voi_ev.size(); ++i)
		std::cout << "VOI " << i + 1 << " dxmc: " << voi_ev[i] << ", TG195: " << sim_ev[i] << ", difference [%]: " << (voi_ev[i] - sim_ev[i]) / sim_ev[i] * 100 << "\n";
	std::cout << "Collimation: 80 mm:\n";
	src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(40. / 600.) * 2.);
	dose = transport::run(w, &src, nullptr, false);

	for (std::size_t i = 0; i < 4; ++i) {
		voi_ev[i] = 1000 / total_hist * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m)->double {return m == i + 2 ? d : 0.0; });
	}
	sim_ev = { 3586.59,3537.84,3378.99,2672.21	};

	for (std::size_t i = 0; i < voi_ev.size(); ++i)
		std::cout << "VOI " << i + 1 << " dxmc: " << voi_ev[i] << ", TG195: " << sim_ev[i] << ", difference [%]: " << (voi_ev[i] - sim_ev[i]) / sim_ev[i] * 100 << "\n";


	std::cout << "\n";

	return false;
}

bool testAttenuation()
{
	const double energy = 56.4;
	Material m("Tissue, Soft (ICRP)");
	m.setStandardDensity(1.3);
	std::array<double, 3> spacing = { 1, 1, 1 };
	std::array<std::size_t, 3> dim = { 1, 1, 200 };
	const auto size = std::accumulate(dim.cbegin(), dim.cend(), 1.0, std::multiplies<>());
	auto dens = std::make_shared<std::vector<double>>(size, m.standardDensity());
	auto mat = std::make_shared<std::vector<unsigned char>>(size, static_cast<unsigned char>(0));
	World w;
	w.setDimensions(dim);
	w.setSpacing(spacing);
	w.setDensityArray(dens);
	w.setMaterialIndexArray(mat);
	w.addMaterialToMap(m);
	w.setAttenuationLutMaxEnergy(std::ceil(energy));
	w.validate();

	PencilSource pen;
	pen.setHistoriesPerExposure(1e7);
	pen.setPhotonEnergy(energy);
	pen.setTotalExposures(16);
	std::array<double, 6> cos = { 1,0,0,0,1,0 };
	pen.setDirectionCosines(cos);
	pen.setPosition(0, 0, -400);

	const auto tot_hist = pen.historiesPerExposure() * pen.totalExposures();

	auto kev = transport::run(w, &pen, nullptr, false);
	std::vector<double> att(kev.size());
	for (int i = 0; i < dim[2]; ++i)
		att[i] = std::exp(-(i + 1) * spacing[2]*0.1 * m.standardDensity() * m.getTotalAttenuation(56.4));
	
	const auto att_max = *std::max_element(att.cbegin(), att.cend());
	const auto kev_max = *std::max_element(kev.cbegin(), kev.cend());

	std::transform(att.cbegin(), att.cend(), att.begin(), [=](double a) {return a / att_max; });
	std::transform(kev.cbegin(), kev.cend(), kev.begin(), [=](double e) {return e / kev_max; });

	const double rms = std::sqrt((1.0 / att.size()) * std::transform_reduce(kev.cbegin(), kev.cend(), att.cbegin(), 0.0, std::plus<>(), [](double e, double a)->double {return (e - a) * (e - a); }));

	std::cout << "Test attenuation for pencil beam in 1mm^2 tissue rod: \n";
	std::cout << "Monochromatic beam of " << energy << " kev. \n";
	std::cout << "RMS differense [%] from analytical attenuation: " << rms*100.0 << "\n";
	/*std::cout << "Data: Position (midtpoint) [mm], dxmc rel. dose, analytical rel. attenuation, Difference [%]\n";
	for (int i = 0; i < dim[2]; ++i)
		std::cout << i * spacing[2] + spacing[2] * .5 << ", " << kev[i] << ", " << att[i] << ", " << (kev[i] - att[i]) / att[i]*100.0 << "\n";
	*/
	if (rms < 0.5) {
		std::cout << "SUCCESS\n\n";
		return true;
	}
	std::cout << "FAILURE\n\n";
	return false;
}


int main(int argc, char* argv[])
{
	auto success = testAttenuation();
	success = success && TG195Case4AbsorbedEnergy();
	success = success && TG195Case2AbsorbedEnergy120();
	success = success && TG195Case2AbsorbedEnergyMono();
	if (success)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}