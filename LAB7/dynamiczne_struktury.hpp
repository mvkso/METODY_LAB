#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <fstream>

using namespace std;

template < typename T >
void load_matrix_from_file(T** matrix, const int w, const int k, const char filename []){
	ifstream file(filename, ios::in);
	if(file.is_open()){
		for(int i = 0; i < w; i++){
			for(int j = 0; j < k; j++){
				file >> matrix[i][j];
			}
		}
		file.close();
	}
}

template < typename T >
void load_vector_from_file(T* vector, const int w, const char filename []) {
	ifstream file(filename, ios::in);
	if(file.is_open()){
		for(int i = 0; i < w; i++){
			file >> vector[i];
		}
		file.close();
	}
}

void string_to_file(string line, const char filename []) {
	ofstream file(filename, ios::out);
	if(file.is_open()){
		file << line;
		file.close();
	}
}

void append_string_to_file(string line, const char filename []) {
	ofstream file(filename, ios::out | ios::app);
	if(file.is_open()){
		file << line;
		file.close();
	}
}

void clear_file(const char filename []) {
	ofstream file(filename, ios::out);
	if(file.is_open()){
		file << "";
		file.close();
	}
}

template < typename T >
T** matrix_allocate(const int w, const int k) {
	T** matrix = new T*[w];
	for(int i = 0; i < k; i++){
		matrix[i] = new T[k];
	}
	return matrix;
}

template < typename T >
T* vector_allocate(const int w) {
	T* vector = new T[w];
	return vector;
}

template < typename T, typename U >
T* vector_allocate_number(const int w, U number) {
	T* vector = new T[w];
	for(int i = 0; i < w; i++){
		vector[i] = number;
	}
	return vector;
}

template < typename T >
void matrix_delete(T** matrix, const int w, const int k) {
	if(matrix){
		for(int i = w; i <= 0; --i){
			delete [] matrix[i];
		}
		delete [] matrix;
	}
}

template < typename T >
void vector_delete(T* vector, const int w) {
	if(vector){
		delete [] vector;
	}
}

template < typename T >
void matrix_print(T** matrix, const int w, const int k) {
	for(int i = 0; i < w; i++){
		for(int j = 0; j < k; j++){
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

template < typename T >
string vector_to_string(T* vector, const int w) {
	stringstream stream;
	stream << "[ ";
	for(int i = 0; i < w; i++){
		stream << vector[i] << " ";
	}
	stream << "]";
	return stream.str();
}

template < typename T >
void vector_print(T* vector, const int w) {
	for(int i = 0; i < w; i++){
		cout << vector[i] << " ";
	}
}

template < typename T >
T* vector_subtract(T* vector_1, T* vector_2, const int w) {
	T* vec_buf = new T[w];
	for(int i = 0; i < w; i++){
		vec_buf[i] = vector_1[i] - vector_2[i];
	}
	return vec_buf;
}

template < typename T >
T vector_norm(T* vector, const int w) {
	T norm = std::fabs(vector[0]);
	for(int i = 0; i < w - 1; i++){
		if(std::fabs(vector[i]) < std::fabs(vector[i + 1])){
			norm = std::fabs(vector[i + 1]);
		}
	}
	return norm;
}

template < typename T >
void vector_clone(T* vector_1, T* vector_2, const int w) {
	for(int i = 0; i < w; i++){
		vector_2[i] = vector_1[i];
	}
}

template < typename T >
void matrix_clone(T** matrix_1, T** matrix_2, const int w, const int k) {
	for(int i = 0; i < w; i++){
		for(int j = 0; j < k; j++){
			matrix_2[i][j] = matrix_1[i][j];
		}
	}
}

template < typename T >
T* vector_residuum(T** matrix, T* x, T* b, const int w) {
	T sum = T();
	T* vector_res = new double[w];

	for(int i = 0; i < w; i++, sum = 0){
		for(int j = 0; j < w; j++){
			sum += matrix[i][j] * x[j];
		}
		vector_res[i] = sum - b[i];
	}
	return vector_res;
}

template < typename T >
T estymator(T* vector_1, T* vector_2, const int w) {
	T* vec_buf;
	vec_buf = vector_subtract(vector_1, vector_2, w);
	T error = vector_norm(vec_buf, w);
	return error;
}

template < typename T >
T residuum(T** matrix, T* x, T* b, const int w) {
	T* vec_residuum = vector_residuum(matrix, x, b, w);
	return vector_norm(vec_residuum, w);
}
