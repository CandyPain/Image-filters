#pragma GCC target("avx2")
#include <iostream>
#include <immintrin.h>
#include <fstream>
#include <png.h>
#include <lodepng.h>
#include <stdio.h>
#include <ctime>
using namespace std;
void decode(string filename, vector<unsigned char>& image, unsigned width, unsigned height);
void encode(string filename, std::vector<unsigned char>& image, unsigned width, unsigned height);
unsigned int start_time;
unsigned int end_time;

void insertionSort(int* arr, int n)
{
	int i, key, j;
	for (i = 1; i < n; i++)
	{
		key = arr[i];
		j = i - 1;
		while (j >= 0 && arr[j] > key)
		{
			arr[j + 1] = arr[j];
			j = j - 1;
		}
		arr[j + 1] = key;
	}
}
void insertionSort_vec(unsigned char* arr, int n)
{
	int i, key, j;
	for (i = 1; i < n; i++)
	{
		key = arr[i];
		j = i - 1;
		while (j >= 0 && arr[j] > key)
		{
			arr[j + 1] = arr[j];
			j = j - 1;
		}
		arr[j + 1] = key;
	}
}

void GetRGB(unsigned char* mas, int ost, vector<unsigned char>& image)
{
	int temp = 0;
	for (int i = ost - 1; i < image.size(); i = i + 4)
	{
		mas[temp] = image[i];
		temp++;
	}
}

void SetRGB(unsigned char* mas, int ost, vector<unsigned char>& image) //ost - 1R, 2G, 3B
{
	int temp = 0;
	for (int i = ost - 1; i < image.size(); i = i + 4)
	{
		image[i] = mas[temp];
		temp++;
	}
}


unsigned char* GetMedian(unsigned char* color, unsigned width, unsigned height)
{
	int window[15 * 15];
	unsigned char* answer = new unsigned char[width * height];
	for (int i = 0; i < width * height; ++i)
	{
		answer[i] = color[i];
	}
	int temp = 0;
	for (int i = 1; i < width * height - width * 15 - 16; ++i)
	{
		temp = 0;
		if (width % i == width - 14)
		{
			i += 14;
		}
		else
		{
			for (int j = i; j < i + 15; ++j)
			{
				for (int k = 0; k < 15; ++k)
				{
					window[temp] = color[j + k * width];
					temp++;
				}
			}
		}
		insertionSort(window, 15 * 15);
		if (i - width < 0)
		{
			for (int k = 0; k < 7; ++k)
			{
				for (int j = i; j < i + 15; ++j)
				{
					answer[j + width * k] = window[15 * 15 / 2];
				}
			}
		}
		if (i + width * 16 > width * height)
		{
			for (int k = 0; k < 7; ++k)
			{
				for (int j = i + width * 7; j < i + width * 7 + 1; ++j)
				{
					answer[j + width * k] = window[15 * 15 / 2];
				}
			}
		}
		answer[i + 7 * width + 7] = window[15 * 15 / 2];
	}
	return answer;
}

unsigned char* GetMedian_vec(unsigned char* color, unsigned width, unsigned height)
{
	unsigned char* window = new unsigned char[16 * 16];
	unsigned char* answer = new unsigned char[width * height];
	for (int i = 0; i < width * height; ++i)
	{
		answer[i] = color[i];
	}
	int temp = 0;
	for (int i = 1; i < width * height - width * 16 - 16; ++i)
	{
		temp = 0;
		if (width % i == width - 15)
		{
			i += 15;
		}
		else
		{
			for (int j = 0; j < 16; j++)
			{
				__m128i* ptr = (__m128i*) window;
				__m128i y = _mm_loadu_si128((__m128i*) & color[i + j * width]);
				_mm_storeu_si128(ptr + j, y);
			}
		}
		insertionSort_vec(window, 16 * 16);
		if (i - width < 0)
		{
			for (int k = 0; k < 7; ++k)
			{
				for (int j = i; j < i + 15; ++j)
				{
					answer[j + width * k] = window[15 * 15 / 2];
				}
			}
		}
		if (i + width * 16 > width * height)
		{
			for (int k = 0; k < 7; ++k)
			{
				for (int j = i + width * 7; j < i + width * 7 + 1; ++j)
				{
					answer[j + width * k] = window[15 * 15 / 2];
				}
			}
		}
		answer[i + 7 * width + 7] = window[16 * 16 / 2];
	}
	delete[] window;
	return answer;
}

unsigned char* GetMedian_omp(unsigned char* color, unsigned width, unsigned height)
{
	unsigned char* answer = new unsigned char[width * height];
	for (int i = 0; i < width * height; ++i)
	{
		answer[i] = color[i];
	}
#pragma omp parallel for
	for (int i = 1; i < width * height - width * 15 - 16; ++i)
	{
		int window[15 * 15];
		int temp = 0;
		if (width % i == width - 14)
		{
			continue;
		}
		else
		{
			for (int j = i; j < i + 15; ++j)
			{
				for (int k = 0; k < 15; ++k)
				{
					window[temp] = color[j + k * width];
					temp++;
				}
			}
		}
		insertionSort(window, 15 * 15);
		if (i - width < 0)
		{
			for (int k = 0; k < 7; ++k)
			{
				for (int j = i; j < i + 15; ++j)
				{
					answer[j + width * k] = window[15 * 15 / 2];
				}
			}
		}
		if (i + width * 16 > width * height)
		{
			for (int k = 0; k < 7; ++k)
			{
				for (int j = i + width * 7; j < i + width * 7 + 1; ++j)
				{
					answer[j + width * k] = window[15 * 15 / 2];
				}
			}
		}
		answer[i + 7 * width + 7] = window[15 * 15 / 2];
	}

	return answer;
}

void Median(vector<unsigned char>& image, string& filename, string& NewFile, unsigned width, unsigned height)
{
	decode(filename, image, width, height);
	start_time = clock();
	unsigned char* R = new unsigned char[width * height];
	unsigned char* G = new unsigned char[width * height];
	unsigned char* B = new unsigned char[width * height];
	GetRGB(R, 1, image);
	GetRGB(G, 2, image);
	GetRGB(B, 3, image);
	R = GetMedian(R, width, height);
	G = GetMedian(G, width, height);
	B = GetMedian(B, width, height);
	SetRGB(R, 1, image);
	SetRGB(G, 2, image);
	SetRGB(B, 3, image);
	end_time = clock();
	encode(NewFile, image, width, height);
}


void Median_vec(vector<unsigned char>& image, string& filename, string& NewFile, unsigned width, unsigned height)
{
	decode(filename, image, width, height);
	start_time = clock();
	unsigned char* R = new unsigned char[width * height];
	unsigned char* G = new unsigned char[width * height];
	unsigned char* B = new unsigned char[width * height];
	GetRGB(R, 1, image);
	GetRGB(G, 2, image); //векторизация ухудшит
	GetRGB(B, 3, image);
	R = GetMedian_vec(R, width, height);
	G = GetMedian_vec(G, width, height);
	B = GetMedian_vec(B, width, height);
	SetRGB(R, 1, image);
	SetRGB(G, 2, image);//векторизация ухудшит
	SetRGB(B, 3, image);
	end_time = clock();
	encode(NewFile, image, width, height);
}

void Median_omp(vector<unsigned char>& image, string& filename, string& NewFile, unsigned width, unsigned height)
{
	decode(filename, image, width, height);
	start_time = clock();
	unsigned char* R = new unsigned char[width * height];
	unsigned char* G = new unsigned char[width * height];
	unsigned char* B = new unsigned char[width * height];
	GetRGB(R, 1, image);
	GetRGB(G, 2, image);
	GetRGB(B, 3, image);
	R = GetMedian_omp(R, width, height);
	G = GetMedian_omp(G, width, height);
	B = GetMedian_omp(B, width, height);
	SetRGB(R, 1, image);
	SetRGB(G, 2, image);
	SetRGB(B, 3, image);
	end_time = clock();
	encode(NewFile, image, width, height);
}


void Negativ(vector<unsigned char>& image, string& filename, string& NewFile, unsigned width, unsigned height)
{
	decode(filename, image, width, height);
	start_time = clock();
	for (int i = 0; i < image.size(); ++i)
	{
		if ((i + 1) % 4 != 0)
		{
			image[i] = 255 - image[i];
		}
	}
	end_time = clock();
	encode(NewFile, image, width, height);
}

void Negativ_vec(vector<unsigned char>& image, string& filename, string& NewFile, unsigned width, unsigned height)
{
	decode(filename, image, width, height);
	start_time = clock();
	int i;
	vector<unsigned char> alfa;
	for (i = 3; i < image.size(); i = i + 4)
	{
		alfa.push_back(image[i]);
	}
	i = 0;
	__m256i* ptr = (__m256i*) image.data();
	__m256i x = _mm256_setzero_si256();
	__m256i mask = _mm256_set1_epi32(255);
	for (i = 0; i + 32 < image.size(); i += 32)
	{
		__m256i y = _mm256_loadu_si256((__m256i*) & image[i]);
		x = _mm256_sub_epi32(mask, y);
		_mm256_storeu_si256(ptr + i / 32, x);
	}
	int temp = 0;
	for (int i = 3; i < image.size(); i = i + 4)
	{
		image[i] = alfa[temp];
		temp++;
	}
	end_time = clock();
	encode(NewFile, image, width, height);
}

void Negativ_omp(vector<unsigned char>& image, string& filename, string& NewFile, unsigned width, unsigned height)
{
	decode(filename, image, width, height);
	start_time = clock();
#pragma omp parallel for
	for (int i = 0; i < image.size(); ++i)
	{
		if ((i + 1) % 4 != 0)
		{
			image[i] = 255 - image[i];
		}
	}
	end_time = clock();
	encode(NewFile, image, width, height);
}

void decode(string filename, vector<unsigned char>& image, unsigned width, unsigned height) {
	unsigned error = lodepng::decode(image, width, height, filename);
	if (error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}

void encode(string filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
	std::vector<unsigned char> png;
	unsigned error = lodepng::encode(png, image, width, height);
	if (!error) lodepng::save_file(png, filename);
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}


int main()
{
	setlocale(LC_ALL, "ru");

	cout << "Введите полный путь до файла" << endl;
	string filename;
	getline(cin, filename);
	string NewFile = filename;
	NewFile.insert(NewFile.find("."), "-answer");
	unsigned width;
	unsigned height;
	cout << "Введите высоту и ширину" << endl;
	cin >> width >> height;
	int filter;
	int parallel;
	cout << "Выберите фильтр\n1-Негатив\n2-Медианный фильтр" << endl;
	cin >> filter;
	cout << "Выберите способ\n1-Последовательно\n2-Open Mp\n3-Векторизация" << endl;
	cin >> parallel;
	vector<unsigned char> image;
	if (filter == 1)
	{
		if (parallel == 1)
		{
			Negativ(image, filename, NewFile, width, height);
		}
		if (parallel == 2)
		{
			Negativ_omp(image, filename, NewFile, width, height);
		}
		if (parallel == 3)
		{
			Negativ_vec(image, filename, NewFile, width, height);
		}
	}
	if (filter == 2)
	{
		if (parallel == 1)
		{
			Median(image, filename, NewFile, width, height);
		}
		if (parallel == 2)
		{
			Median_omp(image, filename, NewFile, width, height);
		}
		if (parallel == 3)
		{
			Median_vec(image, filename, NewFile, width, height);
		}
	}
	cout << "Время обработки составило " << end_time - start_time << " мс.";
}