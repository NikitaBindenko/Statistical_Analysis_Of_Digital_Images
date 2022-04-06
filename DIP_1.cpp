#define _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <Windows.h>

using namespace std;

void readBmp(vector<vector<RGBTRIPLE>>& colours, BITMAPFILEHEADER& bfh, BITMAPINFOHEADER& bih, RGBTRIPLE& rgb, const char* filename) {

    FILE* inputFile = fopen(filename, "r + b");
    fread(&bfh, sizeof(bfh), 1, inputFile);
    fread(&bih, sizeof(bih), 1, inputFile);

    size_t padding = 0;
    if ((bih.biWidth * sizeof(rgb)) % 4) padding = 4 - (bih.biWidth * sizeof(rgb)) % 4;
    int count = 0;
    colours.resize(bih.biHeight);
    for (int i = 0; i < bih.biHeight; i++) {
        colours[i].resize(bih.biWidth);
        for (int j = 0; j < bih.biWidth; j++) {

            fread(&rgb, sizeof(rgb), 1, inputFile);
            colours[i][j] = rgb;
        }
        if (padding != 0) {
            fread(&rgb, padding, 1, inputFile);
        }

    }
    fclose(inputFile);
}

void write(vector<int>& freq, const string& filename, bool isShift = 0, int shift = 0) {
    ofstream out(filename);
    if (!out.is_open())
        throw logic_error("ERROR");
    else {
        for (auto i = 0; i < freq.size(); i++)
            if (isShift) out << i - shift << ' ' << freq[i] << endl;
            else out << i << ' ' << freq[i] << endl;

    }
}

double MathExpectation(vector<vector<double>> A, uint32_t H, uint32_t W) {
    double sum = 0;
    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
            sum += A[i][j];

    return (double)sum / ((double)H * W);
}

double entrophi(vector<int>& x, int n) {
    double H = 0;
    for (int i = 0; i < x.size(); ++i) {
        if (x[i] > 0) {
            H += (double)x[i] / n * log2((double)x[i] / n);
        }
    }
    return -H;
}

double sigma(vector<vector<double>> A, double MathExpectation_A, uint32_t H, uint32_t W) {
    double sum = 0;
    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
            sum += pow(A[i][j] - MathExpectation_A, 2);

    double tmp = (double)1 / ((double)H * W - 1);
    return sqrt(tmp * sum);
}

double CorrelationCoefficient(vector<vector<double>> A, vector<vector<double>> B, uint32_t H, uint32_t W) {
    double MathExpectation_A = MathExpectation(A, H, W);
    double MathExpectation_B = MathExpectation(B, H, W);
    vector<vector<double>> MathExpectation_product(H, vector<double>(W));

    for (int i = 0; i < H; i++)
        for (int j = 0; j < W; j++)
            MathExpectation_product[i][j] = ((A[i][j] - MathExpectation_A) * (B[i][j] - MathExpectation_B));

    double sigma_first = sigma(A, MathExpectation_A, H, W);
    double sigma_second = sigma(B, MathExpectation_B, H, W);
    double result = MathExpectation(MathExpectation_product, H, W) / (sigma_first * sigma_second);
    return result;
}

void colorsToVector(vector<vector<RGBTRIPLE>> colors, vector<vector<double>> &A, int H, int W,  int colorIs) {

    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            if (colorIs == 0) { A[i][j] = colors[i][j].rgbtRed; }
            else if (colorIs == 1) { A[i][j] = colors[i][j].rgbtGreen; }
            else if (colorIs == 2) { A[i][j] = colors[i][j].rgbtBlue; }
            else { break; }
        }
    }
}

vector<vector<double>> FirstSelection(vector<vector<double>> array, int x, int y, int H, int W) {
    vector<vector<double>> result(H - abs(y), vector<double>(W - abs(x)));
    if (y >= 0) {
        if (x >= 0) {
            for (int i = 0; i < array.size() - y; ++i) {
                for (int j = 0; j < array[0].size() - x; ++j) {
                    result[i][j] = array[i][j];
                }
            }
        }
        else {
            for (int i = 0; i < array.size() - y; ++i) {
                for (int j = -x; j < array[0].size(); ++j) {
                    result[i][j + x] = array[i][j];
                }
            }
        }
    }
    else {
        if (x >= 0) {
            for (int i = -y; i < array.size(); ++i) {
                for (int j = 0; j < array[0].size() - x; ++j) {
                    result[i + y][j] = array[i][j];
                }
            }
        }
        else {
            for (int i = -y; i < array.size(); ++i) {
                for (int j = -x; j < array[0].size(); ++j) {
                    result[i + y][j + x] = array[i][j];
                }
            }
        }
    }
    return result;
}

vector<vector<double>> SecondSelection(vector<vector<double>> array, int x, int y, int H, int W) {
    vector<vector<double>> result(H - abs(y), vector<double>(W - abs(x)));;
    if (y >= 0) {
        if (x >= 0) {
            for (int i = y; i < array.size(); ++i) {
                for (int j = x; j < array[0].size(); ++j) {
                    result[i - y][j - x] = array[i][j];
                }
            }
        }
        else {
            for (int i = y; i < array.size(); ++i) {
                for (int j = 0; j < array[0].size() + x; ++j) {
                    result[i - y][j] = array[i][j];
                }
            }
        }
    }
    else {
        if (x >= 0) {
            for (int i = 0; i < array.size() + y; ++i) {
                for (int j = x; j < array[0].size(); ++j) {
                    result[i][j - x] = array[i][j];
                }
            }
        }
        else {
            for (int i = 0; i < array.size() + y; ++i) {
                for (int j = 0; j < array[0].size() + x; ++j) {
                    result[i][j] = array[i][j];
                }
            }
        }
    }
    return result;
}

void autocorrelation(vector<vector<double>> A, string name, int H, int W, int num = 0) {
    ofstream autocorFunc;
    //vector<vector<double>> temp(A.size() - abs(y), vector<double>(A[0].size() - abs(x)));
    //vector<double> autocor(W / 4);
    for (int y = -10; y <= 10; y += 5) {
        autocorFunc.open("autocorrelation_" + name + " " + to_string(y) + ".txt");
        for (int x = num; x < -num + 1; x++) {
            vector<vector<double>> first(H - abs(y), vector<double>(W - abs(x)));;
            vector<vector<double>> second(H - abs(y), vector<double>(W - abs(x)));;
            first = FirstSelection(A, x, y, H, W);
            second = SecondSelection(A, x, y, H, W);
            autocorFunc << x << " " << CorrelationCoefficient(first, second, H - abs(y), W - abs(x));
            autocorFunc << "\n";
        }
        autocorFunc.close();
    }
}

void PSNR(vector<vector<RGBTRIPLE>> colors1, vector<vector<RGBTRIPLE>> colors2, int H, int W) {
    
    double denumeratorRed = 0, denumeratorGreen = 0, denumeratorBlue = 0;
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            denumeratorRed += pow(colors1[i][j].rgbtRed - colors2[i][j].rgbtRed, 2);
            denumeratorGreen += pow(colors1[i][j].rgbtGreen - colors2[i][j].rgbtGreen, 2);
            denumeratorBlue += pow(colors1[i][j].rgbtBlue - colors2[i][j].rgbtBlue, 2);
        }
    }
    double numerator = W * H * pow(pow(2, 8) - 1, 2);

    cout << "PSNR Red:\t" << 10 * log10(numerator / denumeratorRed) << endl;
    cout << "PSNR Green:\t" << 10 * log10(numerator / denumeratorGreen) << endl;
    cout << "PSNR Blue:\t" << 10 * log10(numerator / denumeratorBlue) << endl;
}

void PSNRforComponent(vector<vector<double>> component1, vector<vector<double>> component2, int H, int W, string componentName) {
    double denumerator = 0;
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            denumerator += pow(component1[i][j] - component2[i][j], 2);
        }
    }
    double numerator = W * H * pow(pow(2, 8) - 1, 2);
    cout << componentName << " PSNR:\t" << 10 * log10(numerator / denumerator) << endl;
}

void task3(const char* fileName, const char* outFileName, string onlyColor) {

    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
    RGBTRIPLE rgb;
    vector<vector<RGBTRIPLE>> colors;
    readBmp(colors, bfh, bih, rgb, fileName);

    FILE* outputFile = fopen(outFileName, "w + b");
    fwrite(&bfh, sizeof(bfh), 1, outputFile);
    fwrite(&bih, sizeof(bih), 1, outputFile);

    for (int i = 0; i < bih.biHeight; i++) {
        for (int j = 0; j < bih.biWidth; j++) {
            if (onlyColor != "red") { colors[i][j].rgbtRed = 0; }
            if (onlyColor != "green") { colors[i][j].rgbtGreen = 0; }
            if (onlyColor != "blue") { colors[i][j].rgbtBlue = 0; }
            fwrite(&colors[i][j], sizeof(colors[i][j]), 1, outputFile);
        }
    }
    fclose(outputFile);
}

void task4(const char* fileName) {
    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
    RGBTRIPLE rgb;
    vector<vector<RGBTRIPLE>> colors;
    readBmp(colors, bfh, bih, rgb, fileName);

    double H = bih.biHeight;
    double W = bih.biWidth;

    vector<vector<double>> red(H, vector<double>(W));
    vector<vector<double>> green(H, vector<double>(W));
    vector<vector<double>> blue(H, vector<double>(W));

    colorsToVector(colors, red, H, W, 0);
    colorsToVector(colors, green, H, W, 1);
    colorsToVector(colors, blue, H, W, 2);
    
    double correlation_RG = CorrelationCoefficient(red, green, H, W);
    double correlation_RB = CorrelationCoefficient(red, blue, H, W);
    double correlation_BG = CorrelationCoefficient(green, blue, H, W);

    cout << correlation_RG << endl;
    cout << correlation_RB << endl;
    cout << correlation_BG << endl;

    autocorrelation(red, "red", H, W, -50);
    autocorrelation(green, "green", H, W, -50);
    autocorrelation(blue, "blue", H, W, -50);
}

void task5_and_7(const char* fileName, const char* recoveredFile) {

    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
    RGBTRIPLE rgb;
    vector<vector<RGBTRIPLE>> colors;
    readBmp(colors, bfh, bih, rgb, fileName);
    double H = bih.biHeight;
    double W = bih.biWidth;

    vector<vector<double>> Y(H, vector<double>(W));
    vector<vector<double>> Cb(H, vector<double>(W));
    vector<vector<double>> Cr(H, vector<double>(W));

    for (int i = 0; i < bih.biHeight; ++i) {
        for (int j = 0; j < bih.biWidth; ++j) {
            Y[i][j] = std::round(0.299 * colors[i][j].rgbtRed + 0.587 * colors[i][j].rgbtGreen + 0.114 * colors[i][j].rgbtBlue);
            Cb[i][j] = std::round(0.5643 * (colors[i][j].rgbtBlue - Y[i][j]) + 128);
            Cr[i][j] = std::round(0.7132 * (colors[i][j].rgbtRed - Y[i][j]) + 128);
        }
    }

    double correlation_Y_Cb = CorrelationCoefficient(Y, Cb, H, W);
    double correlation_Y_Cr = CorrelationCoefficient(Y, Cr, H, W);
    double correlation_Cb_Cr = CorrelationCoefficient(Cb, Cr, H, W);
    
    cout << correlation_Y_Cb << endl;
    cout << correlation_Y_Cr << endl;
    cout << correlation_Cb_Cr << endl;

    //autocorrelation(Y, "Y", H, W, -50);
    //autocorrelation(Cb, "Cb", H, W, -50);
    //autocorrelation(Cr, "Cr", H, W, -50);

    FILE* recovered = fopen(recoveredFile, "w + b");      //7 задание
    fwrite(&bfh, sizeof(bfh), 1, recovered);
    fwrite(&bih, sizeof(bih), 1, recovered);

    vector<vector<RGBTRIPLE>> colours;
    colours.resize(bih.biHeight);
    for (int i = 0; i < bih.biHeight; i++) {
        colours[i].resize(bih.biWidth);
    }

    for (int i = 0; i < bih.biHeight; i++) {
        for (int j = 0; j < bih.biWidth; j++) {

            double g = round(Y[i][j] - (0.714 * (Cr[i][j] - 128)) - (0.334 * (Cb[i][j] - 128)));
            double r = round((Y[i][j] + 1.402 * (Cr[i][j] - 128)));
            double b = round((Y[i][j] + 1.772 * (Cb[i][j] - 128)));

                if (g < 0) { g = 0; }   //клиппирование
                if (g > 255) { g = 255; }
                if (r < 0) { r = 0; }
                if (r > 255) { r = 255; }
                if (b < 0) { b = 0; }
                if (b > 255) { b = 255; }

            colours[i][j].rgbtGreen = std::round(g);
            colours[i][j].rgbtRed = std::round(r);
            colours[i][j].rgbtBlue = std::round(b);

            fwrite(&colours[i][j], sizeof(colours[i][j]), 1, recovered);
        }
    }
    PSNR(colors, colours, bih.biHeight, bih.biWidth);
}

void task6(const char* fileName, const char* outFileName, int component) { // Y - intensity[0], Cb - intensity[1], Cr - intensity[2]
    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
    RGBTRIPLE rgb;
    vector<vector<RGBTRIPLE>> colors;
    readBmp(colors, bfh, bih, rgb, fileName);

    FILE* outputFile = fopen(outFileName, "w + b");
    fwrite(&bfh, sizeof(bfh), 1, outputFile);
    fwrite(&bih, sizeof(bih), 1, outputFile);

    double intensity[3] = { 0 };

    for (int i = 0; i < bih.biHeight; ++i) {
        for (int j = 0; j < bih.biWidth; ++j) {
            intensity[0] = std::round(0.299 * colors[i][j].rgbtRed + 0.587 * colors[i][j].rgbtGreen + 0.114 * colors[i][j].rgbtBlue);
            intensity[1]  = std::round(0.5643 * (colors[i][j].rgbtBlue - intensity[0]) + 128);
            intensity[2] = std::round(0.7132 * (colors[i][j].rgbtRed - intensity[0]) + 128);

            colors[i][j].rgbtRed = intensity[component];
            colors[i][j].rgbtGreen = intensity[component];
            colors[i][j].rgbtBlue = intensity[component];
            fwrite(&colors[i][j], sizeof(colors[i][j]), 1, outputFile);
        }
    }
    fclose(outputFile);
}

vector<vector<double>> decimation1(vector<vector<double>> pixel, int H, int W, int coefficient) {
    vector < vector < double>> result(H / coefficient, vector<double>(W / coefficient));
    int x = 0;
    for (int i = 0; i < H; i++) {
        int y = 0;
        for (int j = 0; j < W; j++) {
            if ((j + 1) % coefficient == 0 && (i + 1) % coefficient == 0) {
                result[x][y] = pixel[i][j];
                y++;
                if (y == W / coefficient)x++;
            }
        }
    }
    return result;
}

vector<vector<double>> decimation2(vector<vector<double>> pixel, int H, int W, int coefficient) {
    vector < vector < double>> result(H / coefficient, vector<double>(W / coefficient));
    for (int i = 0; i < H / coefficient; i++) {
        for (int j = 0; j < W / coefficient; j++) {
            result[i][j] = 0;
            for (int k = i * coefficient; k < i * coefficient + coefficient; k++)
                for (int p = j * coefficient; p < j * coefficient + coefficient; p++)
                    result[i][j] += pixel[k][p];
            result[i][j] = round(result[i][j] / ((double)coefficient * coefficient));
        }
    }
    return result;
}

vector<vector<double>> decimation_recovery(vector<vector<double>> pixel, int H, int W, int coefficient) {
    vector<vector<double>> result(H, vector<double>(W));
    for (int i = 0; i < H / coefficient; i++)
        for (int j = 0; j < W / coefficient; j++)
            for (int k = i * coefficient; k < i * coefficient + coefficient; k++)
                for (int p = j * coefficient; p < j * coefficient + coefficient; p++)
                    result[k][p] = pixel[i][j];
    return result;
}

void YCbCr_ToBMP(vector<vector<double>> Y, vector<vector<double>> Cb_recovered, vector<vector<double>> Cr_recovered, vector<vector<RGBTRIPLE>> colors, FILE* outputFile, int H, int W) {
    
    vector<vector<RGBTRIPLE>> newColors(H, vector<RGBTRIPLE>(W));

    for (int i = 0; i < H; i++) {    //восстановление
        for (int j = 0; j < W; j++) {

            double g = round(Y[i][j] - (0.714 * (Cr_recovered[i][j] - 128)) - (0.334 * (Cb_recovered[i][j] - 128)));
            double r = round((Y[i][j] + 1.402 * (Cr_recovered[i][j] - 128)));
            double b = round((Y[i][j] + 1.772 * (Cb_recovered[i][j] - 128)));

            if (g < 0) { g = 0; }   //клиппирование
            if (g > 255) { g = 255; }
            if (r < 0) { r = 0; }
            if (r > 255) { r = 255; }
            if (b < 0) { b = 0; }
            if (b > 255) { b = 255; }

            newColors[i][j].rgbtGreen = std::round(g);
            newColors[i][j].rgbtRed = std::round(r);
            newColors[i][j].rgbtBlue = std::round(b);

            fwrite(&newColors[i][j], sizeof(newColors[i][j]), 1, outputFile);
        }
    }
    PSNR(colors, newColors, H, W);
    fclose(outputFile);
}

void task8_9_10(const char* fileName, const char* outFileName_a, const char* outFileName_b, int coefficient) {
    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
    RGBTRIPLE rgb;
    vector<vector<RGBTRIPLE>> colors;
    readBmp(colors, bfh, bih, rgb, fileName);

    FILE* outputFile_a = fopen(outFileName_a, "w + b");
    fwrite(&bfh, sizeof(bfh), 1, outputFile_a);
    fwrite(&bih, sizeof(bih), 1, outputFile_a);
    FILE* outputFile_b = fopen(outFileName_b, "w + b");
    fwrite(&bfh, sizeof(bfh), 1, outputFile_b);
    fwrite(&bih, sizeof(bih), 1, outputFile_b);

    vector<vector<double>> Y(bih.biHeight, vector<double>(bih.biWidth));
    vector<vector<double>> Cb(bih.biHeight, vector<double>(bih.biWidth));
    vector<vector<double>> Cr(bih.biHeight, vector<double>(bih.biWidth));

    for (int i = 0; i < bih.biHeight; ++i) {    //преобразование в YCbCr и децимация
        for (int j = 0; j < bih.biWidth; ++j) {
            Y[i][j] = std::round(0.299 * colors[i][j].rgbtRed + 0.587 * colors[i][j].rgbtGreen + 0.114 * colors[i][j].rgbtBlue);
            Cb[i][j] = std::round(0.5643 * (colors[i][j].rgbtBlue - Y[i][j]) + 128);
            Cr[i][j] = std::round(0.7132 * (colors[i][j].rgbtRed - Y[i][j]) + 128);
        }
    }
    
    vector<vector<double>> Cb_decimated_a = decimation1(Cb, bih.biHeight, bih.biWidth, coefficient);
    vector<vector<double>> Cr_decimated_a = decimation1(Cr, bih.biHeight, bih.biWidth, coefficient);
    vector<vector<double>> Cb_decimated_b = decimation2(Cb, bih.biHeight, bih.biWidth, coefficient);
    vector<vector<double>> Cr_decimated_b = decimation2(Cr, bih.biHeight, bih.biWidth, coefficient);

    vector<vector<double>> Cb_recovered_a = decimation_recovery(Cb_decimated_a, bih.biHeight, bih.biWidth, coefficient);
    vector<vector<double>> Cr_recovered_a = decimation_recovery(Cr_decimated_a, bih.biHeight, bih.biWidth, coefficient);
    vector<vector<double>> Cb_recovered_b = decimation_recovery(Cb_decimated_b, bih.biHeight, bih.biWidth, coefficient);
    vector<vector<double>> Cr_recovered_b = decimation_recovery(Cr_decimated_b, bih.biHeight, bih.biWidth, coefficient);

    cout << "\na:" << endl;
    PSNRforComponent(Cb, Cb_recovered_a, bih.biHeight, bih.biWidth, "Cb");
    PSNRforComponent(Cr, Cr_recovered_a, bih.biHeight, bih.biWidth, "Cr");
    YCbCr_ToBMP(Y, Cb_recovered_a, Cr_recovered_a, colors, outputFile_a, bih.biHeight, bih.biWidth);
    cout << "\nb:" << endl;
    PSNRforComponent(Cb, Cb_recovered_b, bih.biHeight, bih.biWidth, "Cb");
    PSNRforComponent(Cr, Cr_recovered_b, bih.biHeight, bih.biWidth, "Cr");
    YCbCr_ToBMP(Y, Cb_recovered_b, Cr_recovered_b, colors, outputFile_b, bih.biHeight, bih.biWidth);

}

void task12_13(const char* fileName) {
    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
    RGBTRIPLE rgb;
    vector<vector<RGBTRIPLE>> colors;
    readBmp(colors, bfh, bih, rgb, fileName);

    vector<int>R(pow(2, 8));
    vector<int>G(pow(2, 8));
    vector<int>B(pow(2, 8));
    vector<int>Y(pow(2, 8));
    vector<int>Cb(pow(2, 8));
    vector<int>Cr(pow(2, 8));

    vector<vector<double>> YPixels(bih.biHeight, vector<double>(bih.biWidth));
    vector<vector<double>> CbPixels(bih.biHeight, vector<double>(bih.biWidth));
    vector<vector<double>> CrPixels(bih.biHeight, vector<double>(bih.biWidth));

    for (int i = 0; i < bih.biHeight; ++i) {
        for (int j = 0; j < bih.biWidth; ++j) {
            YPixels[i][j] = std::round(0.299 * colors[i][j].rgbtRed + 0.587 * colors[i][j].rgbtGreen + 0.114 * colors[i][j].rgbtBlue);
            CbPixels[i][j] = std::round(0.5643 * (colors[i][j].rgbtBlue - YPixels[i][j]) + 128);
            CrPixels[i][j] = std::round(0.7132 * (colors[i][j].rgbtRed - YPixels[i][j]) + 128);
        }
    }

    for (int i = 0; i < bih.biHeight; ++i) {
        for (int j = 0; j < bih.biWidth; ++j) {
            ++R[(int)colors[i][j].rgbtRed];
            ++G[(int)colors[i][j].rgbtGreen];
            ++B[(int)colors[i][j].rgbtBlue];
            ++Y[(int)YPixels[i][j]];
            ++Cb[(int)CbPixels[i][j]];
            ++Cr[(int)CrPixels[i][j]];
        }
    }

    double HR = 0, HG = 0, HB = 0, HY = 0, HCb = 0, HCr = 0;
    int n = bih.biHeight * bih.biWidth;
    for (int k = 0; k < 6; ++k) {
        switch (k) {
        case 0:
            for (int i = 0; i < pow(2, 8); ++i) {
                if (R[i] > 0) {
                    HR += (double)R[i] / n * log2((double)R[i] / n);
                }
            }
            HR = -HR;
            break;
        case 1:
            for (int i = 0; i < pow(2, 8); ++i) {
                if (G[i] > 0) {
                    HG += (double)G[i] / n * log2((double)G[i] / n);
                }
            }
            HG = -HG;
            break;
        case 2:
            for (int i = 0; i < pow(2, 8); ++i) {
                if (B[i] > 0) {
                    HB += (double)B[i] / n * log2((double)B[i] / n);
                }
            }
            HB = -HB;
            break;
        case 3:
            for (int i = 0; i < pow(2, 8); ++i) {
                if (Y[i] > 0) {
                    HY += (double)Y[i] / n * log2((double)Y[i] / n);
                }
            }
            HY = -HY;
            break;
        case 4:
            for (int i = 0; i < pow(2, 8); ++i) {
                if (Cb[i] > 0) {
                    HCb += (double)Cb[i] / n * log2((double)Cb[i] / n);
                }
            }
            HCb = -HCb;
            break;
        case 5:
            for (int i = 0; i < pow(2, 8); ++i) {
                if (Cr[i] > 0) {
                    HCr += (double)Cr[i] / n * log2((double)Cr[i] / n);
                }
            }
            HCr = -HCr;
            break;
        }
    }
    cout << "HR = " << HR << " " << "HG = " << HG << " " << "HB = " << HB << " " << "HY = " << HY << " " << "HCb = " << HCb << " " << "HCr = " << HCr << endl;

    write(R, "R.txt");
    write(G, "G.txt");
    write(B, "B.txt");
    write(Y, "Y.txt");
    write(Cb, "Cb.txt");
    write(Cr, "Cr.txt");
}

void task14_15_16(const char* fileName) {

    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
    RGBTRIPLE rgb;
    vector<vector<RGBTRIPLE>> colors;
    readBmp(colors, bfh, bih, rgb, fileName);

    vector<int>R(pow(2, 8));
    vector<int>G(pow(2, 8));
    vector<int>B(pow(2, 8));
    vector<int>Y(pow(2, 8));
    vector<int>Cb(pow(2, 8));
    vector<int>Cr(pow(2, 8));

    vector<vector<double>> YPixels(bih.biHeight, vector<double>(bih.biWidth));
    vector<vector<double>> CbPixels(bih.biHeight, vector<double>(bih.biWidth));
    vector<vector<double>> CrPixels(bih.biHeight, vector<double>(bih.biWidth));

    for (int i = 0; i < bih.biHeight; ++i) {
        for (int j = 0; j < bih.biWidth; ++j) {
            YPixels[i][j] = std::round(0.299 * colors[i][j].rgbtRed + 0.587 * colors[i][j].rgbtGreen + 0.114 * colors[i][j].rgbtBlue);
            CbPixels[i][j] = std::round(0.5643 * (colors[i][j].rgbtBlue - YPixels[i][j]) + 128);
            CrPixels[i][j] = std::round(0.7132 * (colors[i][j].rgbtRed - YPixels[i][j]) + 128);
        }
    }

    vector<int> D1R(2 * pow(2, 8) - 1);
    vector<int> D2R(2 * pow(2, 8) - 1);
    vector<int> D3R(2 * pow(2, 8) - 1);
    vector<int> D4R(2 * pow(2, 8) - 1);

    vector<int> D1G(2 * pow(2, 8) - 1);
    vector<int> D2G(2 * pow(2, 8) - 1);
    vector<int> D3G(2 * pow(2, 8) - 1);
    vector<int> D4G(2 * pow(2, 8) - 1);

    vector<int> D1B(2 * pow(2, 8) - 1);
    vector<int> D2B(2 * pow(2, 8) - 1);
    vector<int> D3B(2 * pow(2, 8) - 1);
    vector<int> D4B(2 * pow(2, 8) - 1);

    vector<int> D1Y(2 * pow(2, 8) - 1);
    vector<int> D2Y(2 * pow(2, 8) - 1);
    vector<int> D3Y(2 * pow(2, 8) - 1);
    vector<int> D4Y(2 * pow(2, 8) - 1);

    vector<int> D1Cb(2 * pow(2, 8) - 1);
    vector<int> D2Cb(2 * pow(2, 8) - 1);
    vector<int> D3Cb(2 * pow(2, 8) - 1);
    vector<int> D4Cb(2 * pow(2, 8) - 1);

    vector<int> D1Cr(2 * pow(2, 8) - 1);
    vector<int> D2Cr(2 * pow(2, 8) - 1);
    vector<int> D3Cr(2 * pow(2, 8) - 1);
    vector<int> D4Cr(2 * pow(2, 8) - 1);

    for (int k = 0; k < 6; ++k) {
        switch (k) {
        case 0:
            for (int l = 0; l < 4; ++l) {
                switch (l) {
                case 0:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = colors[i][j].rgbtRed - colors[i][j - 1].rgbtRed;
                            ++D1R[value + 255];
                        }
                    }
                    break;
                case 1:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = colors[i][j].rgbtRed - colors[i - 1][j].rgbtRed;
                            ++D2R[value + 255];
                        }
                    }
                    break;
                case 2:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = colors[i][j].rgbtRed - colors[i - 1][j - 1].rgbtRed;
                            ++D3R[value + 255];
                        }
                    }
                    break;
                case 3:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int middle = round((colors[i - 1][j].rgbtRed + colors[i - 1][j - 1].rgbtRed + colors[i][j - 1].rgbtRed) / 3);
                            int value = colors[i][j].rgbtRed - middle;
                            ++D4R[value + 255];
                        }
                    }
                    break;
                }
            }
            break;
        case 1:
            for (int l = 0; l < 4; ++l) {
                switch (l) {
                case 0:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = colors[i][j].rgbtGreen - colors[i][j - 1].rgbtGreen;
                            ++D1G[value + 255];
                        }
                    }
                    break;
                case 1:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = colors[i][j].rgbtGreen - colors[i - 1][j].rgbtGreen;
                            ++D2G[value + 255];
                        }
                    }
                    break;
                case 2:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = colors[i][j].rgbtGreen - colors[i - 1][j - 1].rgbtGreen;
                            ++D3G[value + 255];
                        }
                    }
                    break;
                case 3:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int middle = round((colors[i - 1][j].rgbtGreen + colors[i - 1][j - 1].rgbtGreen + colors[i][j - 1].rgbtGreen) / 3);
                            int value = colors[i][j].rgbtGreen - middle;
                            ++D4G[value + 255];
                        }
                    }
                    break;
                }
            }
            break;
        case 2:
            for (int l = 0; l < 4; ++l) {
                switch (l) {
                case 0:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = colors[i][j].rgbtBlue - colors[i][j - 1].rgbtBlue;
                            ++D1B[value + 255];
                        }
                    }
                    break;
                case 1:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = colors[i][j].rgbtBlue - colors[i - 1][j].rgbtBlue;
                            ++D2B[value + 255];
                        }
                    }
                    break;
                case 2:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = colors[i][j].rgbtBlue - colors[i - 1][j - 1].rgbtBlue;
                            ++D3B[value + 255];
                        }
                    }
                    break;
                case 3:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int middle = round((colors[i - 1][j].rgbtBlue + colors[i - 1][j - 1].rgbtBlue + colors[i][j - 1].rgbtBlue) / 3);
                            int value = colors[i][j].rgbtBlue - middle;
                            ++D4B[value + 255];
                        }
                    }
                    break;
                }
            }
            break;
        case 3:
            for (int l = 0; l < 4; ++l) {
                switch (l) {
                case 0:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = YPixels[i][j] - YPixels[i][j - 1];
                            ++D1Y[value + 255];
                        }
                    }
                    break;
                case 1:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = YPixels[i][j] - YPixels[i - 1][j];
                            ++D2Y[value + 255];
                        }
                    }
                    break;
                case 2:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = YPixels[i][j] - YPixels[i - 1][j - 1];
                            ++D3Y[value + 255];
                        }
                    }
                    break;
                case 3:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int middle = round((YPixels[i - 1][j] + YPixels[i - 1][j - 1] + YPixels[i][j - 1]) / 3);
                            int value = YPixels[i][j] - middle;
                            ++D4Y[value + 255];
                        }
                    }
                    break;
                }
            }
            break;
        case 4:
            for (int l = 0; l < 4; ++l) {
                switch (l) {
                case 0:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = CbPixels[i][j] - CbPixels[i][j - 1];
                            ++D1Cb[value + 255];
                        }
                    }
                    break;
                case 1:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = CbPixels[i][j] - CbPixels[i - 1][j];
                            ++D2Cb[value + 255];
                        }
                    }
                    break;
                case 2:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = CbPixels[i][j] - CbPixels[i - 1][j - 1];
                            ++D3Cb[value + 255];
                        }
                    }
                    break;
                case 3:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int middle = round((CbPixels[i - 1][j] + CbPixels[i - 1][j - 1] + CbPixels[i][j - 1]) / 3);
                            int value = CbPixels[i][j] - middle;
                            ++D4Cb[value + 255];
                        }
                    }
                    break;
                }
            }
            break;
        case 5:
            for (int l = 0; l < 4; ++l) {
                switch (l) {
                case 0:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = CrPixels[i][j] - CrPixels[i][j - 1];
                            ++D1Cr[value + 255];
                        }
                    }
                    break;
                case 1:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = CrPixels[i][j] - CrPixels[i - 1][j];
                            ++D2Cr[value + 255];
                        }
                    }
                    break;
                case 2:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int value = CrPixels[i][j] - CrPixels[i - 1][j - 1];
                            ++D3Cr[value + 255];
                        }
                    }
                    break;
                case 3:
                    for (int i = 1; i < bih.biHeight; ++i) {
                        for (int j = 1; j < bih.biWidth; ++j) {
                            int middle = round((CrPixels[i - 1][j] + CrPixels[i - 1][j - 1] + CrPixels[i][j - 1]) / 3);
                            int value = CrPixels[i][j] - middle;
                            ++D4Cr[value + 255];
                        }
                    }
                    break;
                }
            }
            break;
        }
    }

    write(D1R, "D1R.txt", true, 255);
    write(D2R, "D2R.txt", true, 255);
    write(D3R, "D3R.txt", true, 255);
    write(D4R, "D4R.txt", true, 255);

    write(D1G, "D1G.txt", true, 255);
    write(D2G, "D2G.txt", true, 255);
    write(D3G, "D3G.txt", true, 255);
    write(D4G, "D4G.txt", true, 255);

    write(D1B, "D1B.txt", true, 255);
    write(D2B, "D2B.txt", true, 255);
    write(D3B, "D3B.txt", true, 255);
    write(D4B, "D4B.txt", true, 255);

    write(D1Y, "D1Y.txt", true, 255);
    write(D2Y, "D2Y.txt", true, 255);
    write(D3Y, "D3Y.txt", true, 255);
    write(D4Y, "D4Y.txt", true, 255);

    write(D1Cb, "D1Cb.txt", true, 255);
    write(D2Cb, "D2Cb.txt", true, 255);
    write(D3Cb, "D3Cb.txt", true, 255);
    write(D4Cb, "D4Cb.txt", true, 255);

    write(D1Cr, "D1Cr.txt", true, 255);
    write(D2Cr, "D2Cr.txt", true, 255);
    write(D3Cr, "D3Cr.txt", true, 255);
    write(D4Cr, "D4Cr.txt", true, 255);


    double HR = 0, HG = 0, HB = 0, HY = 0, HCb = 0, HCr = 0;
    int n = bih.biHeight * bih.biWidth;
    cout << "Entropy:" << endl;
    HR = entrophi(D1R, n);
    HG = entrophi(D1G, n);
    HB = entrophi(D1B, n);
    HY = entrophi(D1Y, n);
    HCb = entrophi(D1Cb, n);
    HCr = entrophi(D1Cr, n);
    cout << "D1R = " << HR << " " << "D1G = " << HG << " " << "D1B = " << HB << " " << "D1Y = " << HY << " " << "D1Cb = " << HCb << " " << "D1Cr = " << HCr << endl;

    HR = entrophi(D2R, n);
    HG = entrophi(D2G, n);
    HB = entrophi(D2B, n);
    HY = entrophi(D2Y, n);
    HCb = entrophi(D2Cb, n);
    HCr = entrophi(D2Cr, n);
    cout << "D2R = " << HR << " " << "D2G = " << HG << " " << "D2B = " << HB << " " << "D2Y = " << HY << " " << "D2Cb = " << HCb << " " << "D2Cr = " << HCr << endl;

    HR = entrophi(D3R, n);
    HG = entrophi(D3G, n);
    HB = entrophi(D3B, n);
    HY = entrophi(D3Y, n);
    HCb = entrophi(D3Cb, n);
    HCr = entrophi(D3Cr, n);
    cout << "D3R = " << HR << " " << "D3G = " << HG << " " << "D3B = " << HB << " " << "D3Y = " << HY << " " << "D3Cb = " << HCb << " " << "D3Cr = " << HCr << endl;

    HR = entrophi(D4R, n);
    HG = entrophi(D4G, n);
    HB = entrophi(D4B, n);
    HY = entrophi(D4Y, n);
    HCb = entrophi(D4Cb, n);
    HCr = entrophi(D4Cr, n);
    cout << "D4R = " << HR << " " << "D4G = " << HG << " " << "D4B = " << HB << " " << "D4Y = " << HY << " " << "D4Cb = " << HCb << " " << "D4Cr = " << HCr << endl;
}

void extraTask(const char* fileName, const char* outFileName) { // доп
    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
    RGBTRIPLE rgb;
    vector<vector<RGBTRIPLE>> colors;
    readBmp(colors, bfh, bih, rgb, fileName);

    BITMAPINFOHEADER bih_rotated = bih;
    bih_rotated.biWidth = bih.biHeight;
    bih_rotated.biHeight = bih.biWidth;

    FILE* outputFile = fopen(outFileName, "w + b");
    fwrite(&bfh, sizeof(bfh), 1, outputFile);
    fwrite(&bih_rotated, sizeof(bih_rotated), 1, outputFile);

    vector<vector<RGBTRIPLE>> colors_rotated;

    colors_rotated.resize(bih.biWidth);
    for (int i = 0; i < bih.biWidth; i++) {
        colors_rotated[i].resize(bih.biHeight);
    }

    for (int i = 0; i < bih_rotated.biHeight; ++i) {
        for (int j = 0; j < bih_rotated.biWidth; ++j) {

            colors_rotated[i][j] = colors[j][i];
            fwrite(&colors_rotated[i][j], sizeof(colors_rotated[i][j]), 1, outputFile);
        }
    }
    fclose(outputFile);
}

int main() {

    task3("kodim08.bmp", "task3_1.bmp", "red");
    task3("kodim08.bmp", "task3_2.bmp", "green");
    task3("kodim08.bmp", "task3_3.bmp", "blue");

    task4("kodim08.bmp");

    task5_and_7("kodim08.bmp", "recovered.bmp");

    task6("kodim08.bmp", "Y.bmp", 0);
    task6("kodim08.bmp", "Cb.bmp", 1);
    task6("kodim08.bmp", "Cr.bmp", 2);

    cout << "\nFor 2";
    task8_9_10("kodim08.bmp", "decimated_a_2.bmp", "decimated_b_2.bmp", 2);
    cout << "\nFor 4";
    task8_9_10("kodim08.bmp", "decimated_a_4.bmp", "decimated_b_4.bmp", 4);
    

    cout << "tasks 12 - 13" << endl;
    task12_13("kodim08.bmp");

    cout << "tasks 14 - 16" << endl;
    task14_15_16("kodim08.bmp");

    extraTask("kodim08.bmp", "extraTask_270_rotated.bmp");

    return 0;
}