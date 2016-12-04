#include "Twin.h"

std::vector<double> fromFile(const std::string &, Twin &);
std::vector<double> conv2(const std::vector<double> &, const std::vector<double> &);
std::vector<double> zeroInput(Twin &t);

int main() {

	Twin t("Total Response");

	std::vector<double> h = fromFile("Impulse.txt", t);
	std::vector<double> f = fromFile("Input.txt", t);
	std::vector<double> zs = conv2(f, h);
	std::vector<double> zi = zeroInput(t);

	std::vector<double> total;

	for (int i = 0; i < zi.size(); i++)
		total.push_back(zi[i] + zs[i]);
	t.printmulti(total);
	t.toFile();

	t.exit();
	return EXIT_SUCCESS;
}

std::vector<double> fromFile(const std::string &name, Twin &t) {
	double temp;
	std::vector<double> h;

	std::ifstream fin(name);
	if (fin) 
		while (fin >> temp)
			h.push_back(temp);
	else {
		t.println("Unable to open file " + name + ".");
		t.getInput();
	}
	fin.close();

	return h;
}

std::vector<double> conv2(const std::vector<double> &f, const std::vector<double> &h) {
	std::vector<double> y(f.size() + h.size() - 1);
	double sc = .001;

	int i, j, k;

	for (i = 0; i < h.size() - 1; ++i) {
		y[i] = 0;
		for (j = i, k = 0; j >= 0; --j, ++k)
			y[i] += f[j] * h[k];
		y[i] *= sc;
	}

	for (i = h.size() - 1; i < f.size(); ++i) {
		y[i] = 0;
		for (j = i, k = 0; k < h.size(); --j, ++k)
			y[i] += f[j] * h[k];
		y[i] *= sc;
	}

	return y;
}

std::vector<double> zeroInput(Twin &t) {

	mat A(3), I(3), c(3);

	A[0] = { 0, 1, 0 };
	A[1] = { 0, 0, 1 };
	A[2] = { 0, 0, 0 };

	I[0] = { 1, 0, 0 };
	I[1] = { 0, 1, 0 };
	I[2] = { 0, 0, 1 };

	c[0] = { 0, 0, 0 };
	c[1] = { 0, 0, 0 };
	c[2] = { 0, 0, 0 };

	vec a;

	std::vector<std::vector<double>> results;

	double lt, ut, dt;

	while (1) {
		t.println("For (D^3 + a2*D^2 + a1*D + a0) * y(t) = 0, enter a2, a1, and a0.");
		a = t.gvec();
		if (a.size() == 3) break;
		else t.println("Invalid vector.");
	}

	t.println("Enter lower bound of t.");
	t.getInput(lt);

	t.println("Enter upper bound of t.");
	t.getInput(ut);

	while (1) {
		t.println("Enter the initial condition, vector x(");
		t.print(lt);
		t.print(").");
		results.push_back(t.gvec());
		if (a.size() == 3) break;
		else t.println("Invalid vector.");
	}

	t.println("Enter step size of t.");
	t.getInput(dt);

	A[2][0] = a[2] * -1;
	A[2][1] = a[1] * -1;
	A[2][2] = a[0] * -1;

	//
	// Get constant matrix
	//

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			c[i][j] = A[i][j] * dt + I[i][j];

	//
	// Compute result vectors
	//

	for (double r = lt + dt; r <= ut; r += dt) {

		std::vector<double> vtemp;

		for (int i = 0; i < 3; i++) {
			double temp = 0;
			for (int j = 0; j < 3; j++) {
				temp += c[i][j] * results.back()[j];
			}
			vtemp.push_back(temp);
		}

		results.push_back(vtemp);
	}

	//
	// Plot y
	//

	std::vector<double> plot;

	for (auto &i : results)
		plot.push_back(i[0]);

	return plot;

	t.exit();
}