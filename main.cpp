#include<iostream>
#include<fstream>
#include<complex>
#include<cmath>
using namespace std;

#define BYTE unsigned char // 0~255 까지의 unsigned char 자료형을 앞으로 BYTE 라는 단어로 사용 
#define PI 3.141592
#define INPUT_FILE_NAME1 "Lena_gray.bmp"
#define INPUT_FILE_NAME2 "Lena_gray_NOISE.bmp"
#define HEADERSIZE 1078 // Lena_gray.bmp 파일의 헤더의 크기

int main()
{
	// 이미지 저장 코드

	ifstream In_Image; // 파일 읽기
	In_Image.open(INPUT_FILE_NAME1, ios::in | ios::binary); // INPUT_FILE_NAME1 파일을 binary로 읽어옴
	
	// Check if file opened successfully
	if (!In_Image.is_open()) {
		cout << "Error: Cannot open " << INPUT_FILE_NAME1 << endl;
		system("pause");
		return -1;
	}

	int M = 64, N = 64; // 이미지의 크기 정의(단위:픽셀) 
	BYTE* header = new BYTE[HEADERSIZE]; // 이미지의 헤더 정보를 담기 위한 공간 생성
	BYTE** image = new BYTE*[N]; // Lena_gray.bmp 파일의 이미지 데이터를 담기 위한 공간 생성
	BYTE** r = new BYTE*[N]; // 이미지 데이터 중 red에 대한 데이터를 담을 공간 생성
	BYTE** g = new BYTE*[N]; // 이미지 데이터 중 green에 대한 데이터를 담을 공간 생성
	BYTE** b = new BYTE*[N]; // 이미지 데이터 중 blue에 대한 데이터를 담을 공간 생성
	BYTE** result = new BYTE*[N]; // DFT와 IDFT를 통해 노이즈를 제거한 이미지 데이터를 담기 위한 공간 생성
	for (int i = 0; i < N; i++) { // 각각의 데이터를 담을 공간을 2차원의 형태로 만드는 과정(왜냐하면 이미지가 2차원이기 때문)
		image[i] = new BYTE[M * 3]; // 이미지 데이터을 담을 공간의 크기 : N * (M*3) (한 픽셀은 R,G,B 성분으로 이루어져 있으므로 행의 크기는 64*3이 되어야 함) 
		r[i] = new BYTE[M]; // red에 대한 데이터를 담을 공간의 크기 : N * M 
		g[i] = new BYTE[M]; // green에 대한 데이터를 담을 공간의 크기 : N * M 
		b[i] = new BYTE[M]; // blue에 대한 데이터를 담을 공간의 크기 : N * M 
		result[i] = new BYTE[M * 3]; // 노이즈를 제거한 이미지 데이터을 담을 공간의 크기 : N * (M*3) (출력 이미지 생성을 위해 한 픽셀에 R,G,B 성분이 모두 들어가야 함) 
	}

	In_Image.read((char*)header, HEADERSIZE); // 헤더 정보를 HEADERSIZE의 크기만큼 읽어서 header 변수에 저장
	for (int i = 0; i < N; i++) {
		In_Image.read((char*)image[i], 3 * M); // Lena_gray.bmp 파일의 이미지 데이터를 읽고 이를 image 변수에 저장 
	}

	int place;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) { // 이미지 데이터가 저장되어 있는 image 변수에서 R,G,B 값을 복사(실제 BMP 파일에는 B,G,R 순으로 담겨 있음)
			place = 3 * j;
			b[i][j] = image[i][place];  // blue 데이터 
			g[i][j] = image[i][place + 1]; // green 데이터
			r[i][j] = image[i][place + 2];  // red 데이터
		}
	}
	/* 작성해야 하는 코드 부분

	// 2D DFT 함수 코드로 구현
	// 2D DFT 이미지 출력 (Hint : DFT값은 매우 크기 때문에 정규화 작업을 실행해야 합니다)
	// 노이즈 제거 
	// 2D IDFT 함수 코드로 구현

	*/
	/* ========== 2D DFT/IDFT and Noise Removal ========== */

	// Read noise image
	ifstream In_Noise;
	In_Noise.open(INPUT_FILE_NAME2, ios::in | ios::binary);
	
	// Check if noise file opened successfully
	if (!In_Noise.is_open()) {
		cout << "Error: Cannot open " << INPUT_FILE_NAME2 << endl;
		system("pause");
		return -1;
	}
	
	BYTE* noise_header = new BYTE[HEADERSIZE];
	BYTE** noise_image = new BYTE*[N];
	BYTE** noise_r = new BYTE*[N];
	BYTE** noise_g = new BYTE*[N];
	BYTE** noise_b = new BYTE*[N];
	for (int i = 0; i < N; i++) {
		noise_image[i] = new BYTE[M * 3];
		noise_r[i] = new BYTE[M];
		noise_g[i] = new BYTE[M];
		noise_b[i] = new BYTE[M];
	}

	In_Noise.read((char*)noise_header, HEADERSIZE);
	for (int i = 0; i < N; i++) {
		In_Noise.read((char*)noise_image[i], 3 * M);
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			place = 3 * j;
			noise_b[i][j] = noise_image[i][place];
			noise_g[i][j] = noise_image[i][place + 1];
			noise_r[i][j] = noise_image[i][place + 2];
		}
	}

	// Create frequency domain arrays
	complex<double>** F_original = new complex<double>*[N];
	complex<double>** F_noise = new complex<double>*[N];
	for (int i = 0; i < N; i++) {
		F_original[i] = new complex<double>[M];
		F_noise[i] = new complex<double>[M];
	}

	// 2D DFT for original image (using red channel for gray image)
	for (int u = 0; u < N; u++) {
		for (int v = 0; v < M; v++) {
			complex<double> sum(0.0, 0.0);
			for (int x = 0; x < N; x++) {
				for (int y = 0; y < M; y++) {
					double angle = -2.0 * PI * ((double)(u * x) / N + (double)(v * y) / M);
					complex<double> expTerm(cos(angle), sin(angle));
					complex<double> pixel((double)r[x][y], 0.0);
					sum = sum + pixel * expTerm;
				}
			}
			F_original[u][v] = sum;
		}
	}

	// 2D DFT for noise image
	for (int u = 0; u < N; u++) {
		for (int v = 0; v < M; v++) {
			complex<double> sum(0.0, 0.0);
			for (int x = 0; x < N; x++) {
				for (int y = 0; y < M; y++) {
					double angle = -2.0 * PI * ((double)(u * x) / N + (double)(v * y) / M);
					complex<double> expTerm(cos(angle), sin(angle));
					complex<double> pixel((double)noise_r[x][y], 0.0);
					sum = sum + pixel * expTerm;
				}
			}
			F_noise[u][v] = sum;
		}
	}

	// Output DFT magnitude images (normalized)
	BYTE** dft_original = new BYTE*[N];
	BYTE** dft_noise = new BYTE*[N];
	BYTE** dft_original_img = new BYTE*[N];
	BYTE** dft_noise_img = new BYTE*[N];
	for (int i = 0; i < N; i++) {
		dft_original[i] = new BYTE[M];
		dft_noise[i] = new BYTE[M];
		dft_original_img[i] = new BYTE[M * 3];
		dft_noise_img[i] = new BYTE[M * 3];
	}

	// Find max magnitude for normalization
	double max_mag_original = 0.0;
	double max_mag_noise = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			double mag_original = log(1.0 + abs(F_original[i][j]));
			double mag_noise = log(1.0 + abs(F_noise[i][j]));
			if (mag_original > max_mag_original) max_mag_original = mag_original;
			if (mag_noise > max_mag_noise) max_mag_noise = mag_noise;
		}
	}

	// Normalize and shift for visualization
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			int u = (i + N / 2) % N;
			int v = (j + M / 2) % M;
			dft_original[u][v] = (BYTE)(255.0 * log(1.0 + abs(F_original[i][j])) / max_mag_original);
			dft_noise[u][v] = (BYTE)(255.0 * log(1.0 + abs(F_noise[i][j])) / max_mag_noise);
		}
	}

	// Create BMP format for DFT images
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			place = 3 * j;
			dft_original_img[i][place] = dft_original[i][j];
			dft_original_img[i][place + 1] = dft_original[i][j];
			dft_original_img[i][place + 2] = dft_original[i][j];
			dft_noise_img[i][place] = dft_noise[i][j];
			dft_noise_img[i][place + 1] = dft_noise[i][j];
			dft_noise_img[i][place + 2] = dft_noise[i][j];
		}
	}

	// Save DFT images
	ofstream Out_DFT_Original;
	Out_DFT_Original.open("DFT_original.bmp", ios::out | ios::binary);
	if (Out_DFT_Original.is_open()) {
		Out_DFT_Original.write((char*)header, HEADERSIZE);
		for (int i = 0; i < N; i++) {
			Out_DFT_Original.write((char*)dft_original_img[i], 3 * M);
		}
		Out_DFT_Original.close();
		cout << "DFT_original.bmp saved successfully" << endl;
	}

	ofstream Out_DFT_Noise;
	Out_DFT_Noise.open("DFT_noise.bmp", ios::out | ios::binary);
	if (Out_DFT_Noise.is_open()) {
		Out_DFT_Noise.write((char*)noise_header, HEADERSIZE);
		for (int i = 0; i < N; i++) {
			Out_DFT_Noise.write((char*)dft_noise_img[i], 3 * M);
		}
		Out_DFT_Noise.close();
		cout << "DFT_noise.bmp saved successfully" << endl;
	}

	// Noise removal: Apply band-stop filter to remove periodic noise
	// Detect noise peaks in frequency domain and suppress them
	complex<double>** F_filtered = new complex<double>*[N];
	for (int i = 0; i < N; i++) {
		F_filtered[i] = new complex<double>[M];
	}

	// Copy noise DFT to filtered array
	for (int u = 0; u < N; u++) {
		for (int v = 0; v < M; v++) {
			F_filtered[u][v] = F_noise[u][v];
		}
	}

	// Find and suppress noise peaks (high magnitude points excluding DC component)
	double threshold = max_mag_noise / log(1.0 + abs(F_noise[0][0])) * abs(F_noise[0][0]) * 0.3;
	for (int u = 0; u < N; u++) {
		for (int v = 0; v < M; v++) {
			double mag = abs(F_noise[u][v]);
			// Suppress high frequency noise peaks (preserve low frequencies for image content)
			if (mag > threshold && (u > N/8 || v > M/8) && (u < 7*N/8 || v < 7*M/8)) {
				double ratio = 0.1; // Reduce by 90%
				F_filtered[u][v] = F_filtered[u][v] * ratio;
			}
		}
	}

	// 2D IDFT to reconstruct the image
	double** reconstructed = new double*[N];
	for (int i = 0; i < N; i++) {
		reconstructed[i] = new double[M];
	}

	for (int x = 0; x < N; x++) {
		for (int y = 0; y < M; y++) {
			complex<double> sum(0.0, 0.0);
			for (int u = 0; u < N; u++) {
				for (int v = 0; v < M; v++) {
					double angle = 2.0 * PI * ((double)(u * x) / N + (double)(v * y) / M);
					complex<double> expTerm(cos(angle), sin(angle));
					sum = sum + F_filtered[u][v] * expTerm;
				}
			}
			reconstructed[x][y] = real(sum) / (N * M);
		}
	}

	// Clip and convert to BYTE
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			double val = reconstructed[i][j];
			if (val < 0) val = 0;
			if (val > 255) val = 255;
			r[i][j] = (BYTE)val;
			g[i][j] = (BYTE)val;
			b[i][j] = (BYTE)val;
		}
	}

	// Clean up
	for (int i = 0; i < N; i++) {
		delete[] noise_image[i];
		delete[] noise_r[i];
		delete[] noise_g[i];
		delete[] noise_b[i];
		delete[] F_original[i];
		delete[] F_noise[i];
		delete[] F_filtered[i];
		delete[] dft_original[i];
		delete[] dft_noise[i];
		delete[] dft_original_img[i];
		delete[] dft_noise_img[i];
		delete[] reconstructed[i];
	}
	delete[] noise_image;
	delete[] noise_r;
	delete[] noise_g;
	delete[] noise_b;
	delete[] noise_header;
	delete[] F_original;
	delete[] F_noise;
	delete[] F_filtered;
	delete[] dft_original;
	delete[] dft_noise;
	delete[] dft_original_img;
	delete[] dft_noise_img;
	delete[] reconstructed;
	In_Noise.close();

	/* ========== End of DFT/IDFT and Noise Removal ========== */

	// 이미지 생성 코드

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) { // DFT와 IDFT를 통해 노이즈를 제거한 데이터를 result 변수에 하나로 묶기(BMP 파일 생성을 위해 한 픽셀은 24 bit(R,G,B)로 구성되야 함)
			place = 3 * j;
			result[i][place] = b[i][j];  // gray 이미지는 R,G,B 값이 모두 같으므로 하나의 변수만 이용해서 저장해도 됨
			result[i][place + 1] = g[i][j];
			result[i][place + 2] = r[i][j];
		}
	}

	ofstream Out;
	Out.open("result.bmp", ios::out | ios::binary); // 출력 파일 생성
	if (Out.is_open()) {
		Out.write((char*)header, HEADERSIZE); // 출력파일에 이미지의 헤더정보 작성
		for (int i = 0; i < N; i++) {
			Out.write((char*)result[i], 3 * M); // 노이즈를 제거하고 R,G,B이 하나의 픽셀에 묶인 result 데이터를 출력 파일에 작성
		}
		Out.close();
		cout << "result.bmp saved successfully" << endl;
	} else {
		cout << "Error: Cannot create result.bmp" << endl;
	}

	// Memory cleanup
	for (int i = 0; i < N; i++) {
		delete[] image[i];
		delete[] r[i];
		delete[] g[i];
		delete[] b[i];
		delete[] result[i];
	}
	delete[] image;
	delete[] r;
	delete[] g;
	delete[] b;
	delete[] result;
	delete[] header;
	
	In_Image.close();

	system("pause");
	return 0;
}