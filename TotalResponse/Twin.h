#ifndef __T_WIN__
#define __T_WIN__

#include <windows.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

typedef std::vector<double> vec;
typedef std::vector<vec> mat;

class Twin {
public:
	Twin(const std::string &);
	template<typename T>
	void println(const T &);
	void println();
	template<typename T>
	void print(const T &);
	template<typename T>
	void printmulti(const T &);
	template<typename T>
	void getInput(T &);
	void getInput();
	void exit();
	void toFile();
	vec gvec();
	void dvec(const vec&);
	mat gmat();
	void dmat(const mat&);
private:
	void size();
	void line();
	int rows();
	int columns();
	void display();
	void clear();
	void wincheck();
	int x, y;
	std::vector<std::string> buffer;
	std::string pname;

};

void Twin::size() {

	CONSOLE_SCREEN_BUFFER_INFO csbi;

	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
	x = csbi.srWindow.Right - csbi.srWindow.Left + 1;
	y = csbi.srWindow.Bottom - csbi.srWindow.Top + 1;

}

int Twin::rows() {
	size();
	return y;
}

int Twin::columns() {
	size();
	return x;
}

template<typename T>
void Twin::print(const T &in) {

	if (buffer.size() == 0) buffer.push_back("");

	std::ostringstream oss;
	oss << in;

	std::string temp = buffer.back() + oss.str();

	buffer.pop_back();
	buffer.push_back(temp);

}

template<typename T>
void Twin::println(const T &in) {

	std::ostringstream oss;
	oss << in;

	buffer.push_back(oss.str());

}

template<typename T>
void Twin::printmulti(const T &in) {

	for (auto &i : in)
		println(i);

}

void Twin::println() {
	buffer.push_back("");
}

void Twin::display() {

	wincheck();

	line();
	std::cout << "  " << pname << std::endl;
	line();
	std::cout << std::endl;

	for (auto&& i : buffer) {
		if (i != "") std::cout << "  > " << i << std::endl;
		else std::cout << std::endl;
	}

	for (int i = 0; i < y - 7 - int(buffer.size()); i++) std::cout << std::endl;
	line();
	std::cout << "  " << "Input:" << std::endl << "  > ";
	buffer.clear();

}

void Twin::wincheck() {

	size();

	int b = 0;

	for (auto &i : buffer) {

		if (int(i.length()) > b) b = int(i.length());

	}

	while (y - 8 - int(buffer.size()) < 1 || int(pname.length()) > x - 4 || int(b) > x - 8) {

		std::cout << std::endl << "Window too small.";
		std::cout << std::endl << "Enter 'x' to output to file, or resize window and enter 'r': ";

		char temp;
		std::cin >> temp;
		std::cin.ignore(1000, '\n');
		std::cin.clear();

		if (temp == 'x' || temp == 'X') {
			toFile();
			buffer.push_back("[See " + pname + ".txt for output.]");
		}

		size();

	}

}

Twin::Twin(const std::string &name) {

	pname = name;

}

void Twin::line() {

	size();

	std::cout << " ";
	for (int i = 0; i < x - 2; i++) std::cout << "=";
	std::cout << " ";

}

void Twin::clear() {
	buffer.clear();
}

template<typename T>
void Twin::getInput(T &var) {
	std::vector<std::string> backup = buffer;

	while (true) {
		display();

		if (std::cin >> var) {
			std::cin.ignore(1000, '\n');
			std::cin.clear();

			buffer.push_back("Input:");
			buffer.push_back("");
			println(var);
			buffer.push_back("");
			buffer.push_back("Enter 'x' to save and continue.");
			display();

			char temp;
			std::cin >> temp;
			std::cin.ignore(1000, '\n');
			std::cin.clear();

			if (temp == 'x' || temp == 'X') break;
			else buffer = backup;
		}
		else {

			std::cin.clear();
			std::cin.ignore(1000, '\n');

			buffer.push_back("Invalid entry.");
			buffer.push_back("");
			buffer.push_back("Press enter to continue.");
			display();
			buffer = backup;

			std::cin.ignore(1000, '\n');
			std::cin.clear();
		}
	}

}

void Twin::getInput() {

	buffer.push_back("Press enter to continue.");
	display();

	std::cin.ignore(1000, '\n');
	std::cin.clear();

}


void Twin::exit() {

	buffer.push_back("Press enter to exit. . .");
	display();

	std::cin.ignore(1000, '\n');
	std::cin.clear();

}

void Twin::toFile() {

	std::ofstream fout(pname + ".txt");
	if (fout.is_open()) {
		for (auto &i : buffer)
			fout << i << std::endl;

		fout.close();
		buffer.clear();
		buffer.push_back("Saved output to " + pname + ".txt.");
		buffer.push_back("");
		getInput();
	}
	else std::cout << std::endl << "Unable to open file.";
}

vec Twin::gvec() {

	display();
	vec b;
	double in;
	std::string input;

	std::getline(std::cin, input);
	std::stringstream iss(input);
	iss >> in;
	b.push_back(in);

	while (true) {

		println("Vector:");
		println();
		dvec(b);
		println("Enter the next element.");
		println("Enter 'x' to save and continue.");
		display();

		std::getline(std::cin, input);
		std::stringstream iss(input);
		iss >> in;

		if (input == "x" || input == "X") return b;
		else b.push_back(in);

	}
}

void Twin::dvec(const vec& vec) {

	for (auto&& i : vec) {
		println("[ ");
		print(i);
		print(" ]");
	}

	println();
}

mat Twin::gmat() {

	display();
	mat a;
	double in;
	std::string input;

	std::getline(std::cin, input);
	std::stringstream iss(input);
	a.push_back(vec());
	while (iss >> in) a[0].push_back(in);

	for (int i = 1; i > 0; i++) {

		println("Matrix:");
		println();
		dmat(a);
		println("Enter the next row.");
		println("Enter 'x' to save and continue.");
		display();

		std::getline(std::cin, input);
		std::stringstream iss(input);
		if (input == "x" || input == "X") i = -1;

		else {
			a.push_back(vec());
			while (iss >> in) a[i].push_back(in);

			if (a[i].size() != a[0].size()) {
				println("Invalid row.");
				a.pop_back();
				i--;
			}
		}
	}
	return a;
}

void Twin::dmat(const mat& mat) {

	println();

	for (int i = 0; i < int(mat.size()); i++) {

		print("[ ");

		for (int j = 0; j < int(mat[0].size()); j++) {
			print(mat[i][j]);
			if (j < int(mat[0].size()) - 1) print(",");
			print(" ");
		}

		print("]");
		println();

	}
}
#endif