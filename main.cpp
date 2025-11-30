#include<iostream>
#include<fstream>
#include<complex>
#include<cmath>
using namespace std;

#define BYTE unsigned char // 0~255 ������ unsigned char �ڷ����� ������ BYTE ��� �ܾ�� ��� 
#define PI 3.141592
#define INPUT_FILE_NAME1 "Lena_gray.bmp"
#define INPUT_FILE_NAME2 "Lena_gray_NOISE.bmp"
#define HEADERSIZE 1078 // Lena_gray.bmp ������ ����� ũ��

void main()
{
	// �̹��� ���� �ڵ�

	ifstream In_Image; // ���� �б�
	In_Image.open(INPUT_FILE_NAME1, ios::in | ios::binary); // INPUT_FILE_NAME1 ������ binary�� �о��

	int M = 64, N = 64; // �̹����� ũ�� ����(����:�ȼ�) 
	BYTE* header = new BYTE[HEADERSIZE]; // �̹����� ��� ������ ��� ���� ���� ����
	BYTE** image = new BYTE * [N]; // Lena_gray.bmp ������ �̹��� �����͸� ��� ���� ���� ����
	BYTE** r = new BYTE * [N]; // �̹��� ������ �� red�� ���� �����͸� ���� ���� ����
	BYTE** g = new BYTE * [N]; // �̹��� ������ �� green�� ���� �����͸� ���� ���� ����
	BYTE** b = new BYTE * [N]; // �̹��� ������ �� blue�� ���� �����͸� ���� ���� ����
	BYTE** result = new BYTE * [N]; // DFT�� IDFT�� ���� ����� ������ �̹��� �����͸� ��� ���� ���� ����
	for (int i = 0; i < N; i++) { // ������ �����͸� ���� ������ 2������ ���·� ����� ����(�ֳ��ϸ� �̹����� 2�����̱� ����)
		image[i] = new BYTE[M * 3]; // �̹��� �������� ���� ������ ũ�� : N * (M*3) (�� �ȼ��� R,G,B �������� �̷���� �����Ƿ� ���� ũ��� 64*3�� �Ǿ�� ��) 
		r[i] = new BYTE[M]; // red�� ���� �����͸� ���� ������ ũ�� : N * M 
		g[i] = new BYTE[M]; // green�� ���� �����͸� ���� ������ ũ�� : N * M 
		b[i] = new BYTE[M]; // blue�� ���� �����͸� ���� ������ ũ�� : N * M 
		result[i] = new BYTE[M * 3]; // ����� ������ �̹��� �������� ���� ������ ũ�� : N * (M*3) (��� �̹��� ������ ���� �� �ȼ��� R,G,B ������ ��� ���� ��) 
	}

	In_Image.read((char*)header, HEADERSIZE); // ��� ������ HEADERSIZE�� ũ�⸸ŭ �о header ������ ����
	for (int i = 0; i < N; i++) {
		In_Image.read((char*)image[i], 3 * M); // Lena_gray.bmp ������ �̹��� �����͸� �а� �̸� image ������ ���� 
	}

	int place;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) { // �̹��� �����Ͱ� ����Ǿ� �ִ� image �������� R,G,B ���� ����(���� BMP ���Ͽ��� B,G,R ������ ��� ����)
			place = 3 * j;
			b[i][j] = image[i][place];  // blue ������ 
			g[i][j] = image[i][place + 1]; // green ������
			r[i][j] = image[i][place + 2];  // red ������
		}
	}
	/* �ۼ��ؾ� �ϴ� �ڵ� �κ�

	// 2D DFT �Լ� �ڵ�� ����
	// 2D DFT �̹��� ��� (Hint : DFT���� �ſ� ũ�� ������ ����ȭ �۾��� �����ؾ� �մϴ�)
	// ������ ����
	// 2D IDFT �Լ� �ڵ�� ����

	*/
	/* ========== 2D DFT/IDFT and Noise Removal ========== */

	cout << "========== 2D DFT/IDFT and Noise Removal Process Started ==========" << endl;

	// Read noise image
	cout << "1. Reading noise image..." << endl;
	ifstream In_Noise;
	In_Noise.open(INPUT_FILE_NAME2, ios::in | ios::binary);
	BYTE* noise_header = new BYTE[HEADERSIZE];
	BYTE** noise_image = new BYTE * [N];
	BYTE** noise_r = new BYTE * [N];
	BYTE** noise_g = new BYTE * [N];
	BYTE** noise_b = new BYTE * [N];
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
	cout << "   - Noise image loaded successfully." << endl;

	// Create frequency domain arrays
	complex<double>** F_original = new complex<double>*[N];
	complex<double>** F_noise = new complex<double>*[N];
	for (int i = 0; i < N; i++) {
		F_original[i] = new complex<double>[M];
		F_noise[i] = new complex<double>[M];
	}

	// 2D DFT for original image (using red channel for gray image)
	cout << "2. Performing 2D DFT on original image..." << endl;
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
	cout << "   - Original image DFT completed." << endl;

	// 2D DFT for noise image
	cout << "3. Performing 2D DFT on noise image..." << endl;
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
	cout << "   - Noise image DFT completed." << endl;

	// Output DFT magnitude images (normalized)
	cout << "4. Normalizing and saving DFT images..." << endl;
	BYTE** dft_original = new BYTE * [N];
	BYTE** dft_noise = new BYTE * [N];
	BYTE** dft_original_img = new BYTE * [N];
	BYTE** dft_noise_img = new BYTE * [N];
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
	Out_DFT_Original.write((char*)header, HEADERSIZE);
	for (int i = 0; i < N; i++) {
		Out_DFT_Original.write((char*)dft_original_img[i], 3 * M);
	}
	Out_DFT_Original.close();

	ofstream Out_DFT_Noise;
	Out_DFT_Noise.open("DFT_noise.bmp", ios::out | ios::binary);
	Out_DFT_Noise.write((char*)noise_header, HEADERSIZE);
	for (int i = 0; i < N; i++) {
		Out_DFT_Noise.write((char*)dft_noise_img[i], 3 * M);
	}
	Out_DFT_Noise.close();
	cout << "   - DFT images saved (DFT_original.bmp, DFT_noise.bmp)." << endl;

	// Noise removal: Detect and suppress noise peaks in frequency domain
	cout << "5. Applying noise removal filter..." << endl;
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

	// Compare original and noise DFT to find noise-specific peaks
	// Calculate magnitude difference between noise and original
	struct NoisePeak {
		int u, v;
		double mag_noise;
		double mag_original;
		double difference;
		double ratio;
	};
	
	NoisePeak* noise_peaks = new NoisePeak[N * M];
	int peak_count = 0;
	
	for (int u = 0; u < N; u++) {
		for (int v = 0; v < M; v++) {
			// Skip DC component
			if (u == 0 && v == 0) continue;
			
			double mag_orig = abs(F_original[u][v]);
			double mag_noise = abs(F_noise[u][v]);
			
			noise_peaks[peak_count].u = u;
			noise_peaks[peak_count].v = v;
			noise_peaks[peak_count].mag_noise = mag_noise;
			noise_peaks[peak_count].mag_original = mag_orig;
			noise_peaks[peak_count].difference = mag_noise - mag_orig;
			
			// Calculate ratio (how much bigger is noise magnitude)
			if (mag_orig > 1.0) {
				noise_peaks[peak_count].ratio = mag_noise / mag_orig;
			} else {
				noise_peaks[peak_count].ratio = mag_noise / 1.0;
			}
			
			peak_count++;
		}
	}
	
	// Sort by difference (descending) to find biggest noise contributions
	for (int i = 0; i < peak_count - 1; i++) {
		for (int j = 0; j < peak_count - i - 1; j++) {
			if (noise_peaks[j].difference < noise_peaks[j + 1].difference) {
				NoisePeak temp = noise_peaks[j];
				noise_peaks[j] = noise_peaks[j + 1];
				noise_peaks[j + 1] = temp;
			}
		}
	}
	
	// Find threshold based on top differences
	int percentile_idx = (int)(peak_count * 0.001);  // Top 0.1%
	if (percentile_idx < 10) percentile_idx = 10;
	double diff_threshold = noise_peaks[percentile_idx].difference;
	
	cout << "   - Top noise peak: mag_noise=" << noise_peaks[0].mag_noise 
		 << ", mag_orig=" << noise_peaks[0].mag_original 
		 << ", diff=" << noise_peaks[0].difference 
		 << " at (" << noise_peaks[0].u << "," << noise_peaks[0].v << ")" << endl;
	cout << "   - Difference threshold: " << diff_threshold << endl;
	
	// Apply notch filter to peaks where noise magnitude >> original magnitude
	int peaks_removed = 0;
	
	for (int i = 0; i < peak_count; i++) {
		// Stop if difference becomes insignificant
		if (noise_peaks[i].difference < diff_threshold * 0.5) break;
		
		int u = noise_peaks[i].u;
		int v = noise_peaks[i].v;
		double mag_orig = noise_peaks[i].mag_original;
		double mag_noise = noise_peaks[i].mag_noise;
		
		// Only suppress if noise is significantly larger (ratio > 2)
		if (noise_peaks[i].ratio > 2.0 || noise_peaks[i].difference > diff_threshold) {
			// Calculate suppression factor: reduce to original level
			double target_reduction = 0.0;
			if (mag_noise > 0) {
				target_reduction = 1.0 - (mag_orig / mag_noise);
				// Limit suppression to avoid artifacts
				if (target_reduction > 0.95) target_reduction = 0.95;
				if (target_reduction < 0.3) target_reduction = 0.3;
			}
			
			// Apply notch filter around this peak
			int radius = 2;
			for (int du = -radius; du <= radius; du++) {
				for (int dv = -radius; dv <= radius; dv++) {
					int nu = u + du;
					int nv = v + dv;
					if (nu >= 0 && nu < N && nv >= 0 && nv < M) {
						double dist = sqrt((double)(du * du + dv * dv));
						if (dist <= radius) {
							// Gaussian weight: stronger at center
							double weight = exp(-dist * dist / 2.0);
							double suppression = target_reduction * weight;
							F_filtered[nu][nv] = F_filtered[nu][nv] * (1.0 - suppression);
						}
					}
				}
			}
			peaks_removed++;
		}
	}
	
	delete[] noise_peaks;
	cout << "   - Noise-specific peaks detected and suppressed: " << peaks_removed << " peaks." << endl;

	// 2D IDFT to reconstruct the image
	cout << "6. Performing 2D IDFT to reconstruct image..." << endl;
	double** reconstructed = new double* [N];
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
	cout << "   - IDFT completed." << endl;

	// Clip and convert to BYTE
	cout << "7. Converting reconstructed image to output format..." << endl;
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
	cout << "   - Image conversion completed." << endl;

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

	cout << "========== Process Completed Successfully! ==========" << endl;
	cout << "Output files created:" << endl;
	cout << "   - DFT_original.bmp : DFT of original image" << endl;
	cout << "   - DFT_noise.bmp    : DFT of noise image" << endl;
	cout << "   - result.bmp       : Denoised image" << endl;

	/* ========== End of DFT/IDFT and Noise Removal ========== */

	// �̹��� ���� �ڵ�

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) { // DFT�� IDFT�� ���� ����� ������ �����͸� result ������ �ϳ��� ����(BMP ���� ������ ���� �� �ȼ��� 24 bit(R,G,B)�� �����Ǿ� ��)
			place = 3 * j;
			result[i][place] = b[i][j];  // gray �̹����� R,G,B ���� ��� �����Ƿ� �ϳ��� ������ �̿��ؼ� �����ص� ��
			result[i][place + 1] = g[i][j];
			result[i][place + 2] = r[i][j];
		}
	}

	ofstream Out;
	Out.open("result.bmp", ios::out | ios::binary); // ��� ���� ����
	Out.write((char*)header, HEADERSIZE); // ������Ͽ� �̹����� ������� �ۼ�
	for (int i = 0; i < N; i++) {
		Out.write((char*)result[i], 3 * M); // ����� �����ϰ� R,G,B�� �ϳ��� �ȼ��� ���� result �����͸� ��� ���Ͽ� �ۼ�
	}

	system("pause");
	return;
}