#include "align.h"
#include <string>

#define _USE_MATH_DEFINES
#include <cmath> 
#include <queue>

using std::string;
using std::cout;
using std::endl;
using std::tuple;
using std::make_tuple;
using std::tie;
using std::max;
using std::min;
using std::queue;

#define kDRow 15
#define kDCol 15
#define kCutMagic 0

static const int kDx[8] = {0, -1, -1, -1, 0, 1, 1, 1};
static const int kDy[8] = {1, 1, 0, -1, -1, -1, 0, 1};

enum {
    USE_R_CHANNEL = 4, 
    USE_G_CHANNEL = 2,
    USE_B_CHANNEL = 1
};

static bool useMirror = false;

static uint brightness(tuple<uint, uint, uint> color)
{
    uint r, g, b;

    tie(r, g, b) = color;

    return (0.2126 * r + 0.7152 * g + 0.0722 * b);
}

static int G_value(Image &Kx, Image &Ky, int row, int col) 
{
    int n_rows = Kx.n_rows;
    int n_cols = Kx.n_cols;

    if (row < 0 || col < 0 || row >= n_rows || col >= n_cols)
        return 0;

    uint rx, gx, bx;
    uint ry, gy, by;

    tie(rx, gx, bx) = Kx(row, col);
    tie(ry, gy, by) = Ky(row, col);

    return sqrt(bx * bx + by * by);
}

static long double calculateMSE(
    Image img1, Image img2, 
    int row1, int col1,
    int row2, int col2,
    int n_rows, int n_cols) 
{
    long long sum = 0;
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            int c1, c2;

            c1 = brightness(img1(row1 + i, col1 + j));
            c2 = brightness(img2(row2 + i, col2 + j));

            sum += (c1 - c2) * (c1 - c2);
        }
    }

    long double result = sum / (1.0 * n_rows * n_cols);

    return result;
}

static long double calculateCC(
    Image img1, Image img2, 
    int row1, int col1,
    int row2, int col2,
    int n_rows, int n_cols) 
{
    long long sum = 0;
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            int c1, c2;

            c1 = brightness(img1(row1 + i, col1 + j));
            c2 = brightness(img2(row2 + i, col2 + j));

            sum += c1 * c2;
        }
    }

    long double result = sum;

    return result;
}

static void align_two_image(
    Image &img1, Image &img2, 
    int mask1, int mask2, 
    Image &outImage, bool useMSE) 
{
    if (img1.n_rows != img2.n_rows || img1.n_cols != img2.n_cols)
        return;

    int n_rows  = img1.n_rows;
    int n_cols = img1.n_cols;

    long double metric = -1.0;
    int m_row1 = 0, m_col1 = 0, m_row2 = 0, m_col2 = 0, m_n_rows = 0, m_n_cols = 0;
    for (int dRow = -kDRow; dRow <= kDRow; dRow++) {
        for (int dCol = -kDCol; dCol <= kDCol; dCol++) {
            int c_row1 = max(0, dRow); 
            int c_col1 = max(0, dCol);
            int c_row2 = -min(0, dRow);
            int c_col2 = -min(0, dCol);
            int c_n_rows = n_rows - abs(dRow);
            int c_n_cols = n_cols - abs(dCol);

            long double c_metric;

            if (useMSE) {
                c_metric = calculateMSE(
                    img1, img2, 
                    c_row1, c_col1, c_row2, c_col2, c_n_rows, c_n_cols);
            } else {
                c_metric = calculateCC(
                    img1, img2, 
                    c_row1, c_col1, c_row2, c_col2, c_n_rows, c_n_cols);
            }

            if (metric < 0 || 
                (useMSE  && c_metric < metric) ||
                (!useMSE && c_metric > metric)) {
                metric = c_metric;
                m_row1 = c_row1;
                m_col1 = c_col1;
                m_row2 = c_row2;
                m_col2 = c_col2;
                m_n_rows = c_n_rows;
                m_n_cols = c_n_cols;
            }
        }
    }

    for (int i = 0; i < m_n_rows; i++) {
        for (int j = 0; j < m_n_cols; j++) {
            uint r1, g1, b1, r2, g2, b2;

            tie(r1, g1, b1) = img1(i + m_row1, j + m_col1);
            tie(r2, g2, b2) = img2(i + m_row2, j + m_col2);

            uint r, g, b;
            r = g = b = 0;

            if (mask1 & USE_R_CHANNEL)
                r = r1;
            if (mask1 & USE_G_CHANNEL)
                g = g1;
            if (mask1 & USE_B_CHANNEL) {
                b = b1;
            }

            if (mask2 & USE_R_CHANNEL)
                r = r2;
            if (mask2 & USE_G_CHANNEL)
                g = g2;
            if (mask2 & USE_B_CHANNEL)
                b = b2;

            outImage(i, j) = make_tuple(r, g, b);
        }
    }
}

// TODO:
Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{
    useMirror = isMirror;

    Matrix<double> kernel = {{0, 0, 0},
                             {0, 1, 0},
                             {0, 0, 0}};

    Image extImage = custom(srcImage, kernel);

    return extImage;

    ///*
    Image cannyImage = canny(srcImage, 5, 75);

    int x, y, rows, cols;

    cut_frame(cannyImage, x, y, rows, cols);

    Image cropImage(rows, cols);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cropImage(i, j) = srcImage(i + x, j + y);
        }
    }

    save_image(cropImage, "/home/alexmihy/mash1718/task1/template/out/crop.bmp");

    srcImage = cropImage;
    //*/
    int n_rows  = srcImage.n_rows;
    int n_cols = srcImage.n_cols;
    int out_n_rows = n_rows / 3;
    int out_n_cols = n_cols;

    Image rImage(out_n_rows, out_n_cols);
    Image gImage(out_n_rows, out_n_cols);
    Image bImage(out_n_rows, out_n_cols);

    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < out_n_rows; i++) {
            for (int j = 0; j < out_n_cols; j++) {
                uint br;

                br = brightness(srcImage(k * out_n_rows + i, j));

                switch (k) {
                case 0:
                    rImage(i, j) = make_tuple(br, br, br);
                    break;
                case 1:
                    gImage(i, j) = make_tuple(br, br, br);
                    break;
                case 2:
                    bImage(i, j) = make_tuple(br, br, br);
                    break;
                default:
                    break;
                }
            }
        }
    }

    save_image(rImage, "/home/alexmihy/mash1718/task1/template/out/R.bmp");
    save_image(gImage, "/home/alexmihy/mash1718/task1/template/out/G.bmp");
    save_image(bImage, "/home/alexmihy/mash1718/task1/template/out/B.bmp");

    Image rgImage(out_n_rows, out_n_cols);
    Image outImage(out_n_rows, out_n_cols);

    align_two_image(rImage, gImage, USE_R_CHANNEL, USE_G_CHANNEL, rgImage, true);
    save_image(rgImage, "/home/alexmihy/mash1718/task1/template/out/RG.bmp");

    align_two_image(rgImage, bImage, USE_R_CHANNEL | USE_G_CHANNEL, USE_B_CHANNEL, outImage, true);

    return outImage;
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
    Matrix<double> kernel = {{-1.0 / 6, -2.0 / 3, -1.0 / 6}, 
                             {-2.0 / 3, 13.0 / 3, -2.0 / 3}, 
                             {-1.0 / 6, -2.0 / 3, -1.0 / 6}};

    return custom(src_image, kernel);
}

Image gray_world(Image src_image) {
    int n_rows = src_image.n_rows;
    int n_cols = src_image.n_cols;

    long long r_sum, g_sum, b_sum;
    r_sum = g_sum = b_sum = 0;

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            uint r, g, b;

            tie(r, g, b) = src_image(i, j);

            r_sum += r;
            g_sum += g;
            b_sum += b;
        }
    }

    int r_avg, g_avg, b_avg, avg;

    r_avg = r_sum / (n_rows * n_cols);
    g_avg = g_sum / (n_rows * n_cols);
    b_avg = b_sum / (n_rows * n_cols);
    avg = (r_avg + g_avg + b_avg) / 3;

    Image outImage(n_rows, n_cols);

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            int r, g, b;

            tie(r, g, b) = src_image(i, j);

            outImage(i, j) = make_tuple(
                max(min(r * avg / r_avg, 255), 0), 
                max(min(g * avg / g_avg, 255), 0), 
                max(min(b * avg / b_avg, 255), 0));

        }
    }

    return outImage;
}

Image resize(Image src_image, double scale) {
    return src_image;
}

Image custom(Image src_image, Matrix<double> kernel) {
    int n_rows = src_image.n_rows;
    int n_cols = src_image.n_cols;
    int rad = kernel.n_rows / 2;

    Image extImage;

    if (useMirror) {
        extImage = Image(n_rows + 2 * rad, n_cols + 2 * rad);

        for (int i = 0; i < n_rows; i++)
            for (int j = 0; j < n_cols; j++)
                extImage(i + rad, j + rad) = src_image(i, j);

        for (int k = 0; k < rad; k++) {
            for (int i = 0; i < n_rows; i++) {
                extImage(i + rad, k) = src_image(i, rad - k);
                extImage(i + rad, n_cols + rad + k) = src_image(i, (n_cols - 1) - k);
            }

            for (int j = 0; j < n_cols; j++) {
                extImage(k, j + rad) = src_image(rad - k, j);
                extImage(n_rows + rad + k, j + rad) = src_image((n_rows - 1) - k, j);
            }
        }

        for (int i = 0; i < rad; i++) {
            for (int j = 0; j < rad; j++) {
                extImage(i, j) = src_image(rad - i, rad - j);
                extImage(i + n_rows + rad, j) = src_image((n_rows - 1) - i, rad - j);
                extImage(i, j + n_cols + rad) = src_image(rad - i, (n_cols - 1) - j);
                extImage(i + n_rows + rad, j + n_cols + rad) = src_image((n_rows - 1) - i, (n_cols - 1) - j);
            }
        }
    } else {
        extImage = src_image.deep_copy();
    }

    n_rows = extImage.n_rows - 2 * rad;
    n_cols = extImage.n_cols - 2 * rad;

    Image outImage(n_rows, n_cols);

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            int R, G, B;
            double r_sum, g_sum, b_sum;

            r_sum = g_sum = b_sum = 0;
            for (int k1 = -rad; k1 <= rad; k1++) {
                for (int k2 = -rad; k2 <= rad; k2++) {
                    uint r, g, b;

                    tie(r, g, b) = extImage(i + rad + k1, j + rad + k2);


                    r_sum += r * kernel(k1 + rad, k2 + rad);
                    g_sum += g * kernel(k1 + rad, k2 + rad);
                    b_sum += b * kernel(k1 + rad, k2 + rad);
                }   
            }

            R = r_sum;
            G = g_sum; 
            B = b_sum;

            outImage(i, j) = make_tuple(
                max(min(R, 255), 0), 
                max(min(G, 255), 0), 
                max(min(B, 255), 0));
        }
    }

    return outImage;
}

Image autocontrast(Image src_image, double fraction) {
    return src_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
    int kernel_size = 2 * radius + 1;
    Matrix<double> kernel(kernel_size, kernel_size);

    long double sum = 0.0;
    for (int i = -radius; i <= radius; i++) {
        for (int j = -radius; j <= radius; j++) {
            double value;

            value = 1.0 / (2 * M_PI * sigma * sigma * exp((i * i + j * j) / (2 * sigma * sigma)));
            sum = sum + value;

            kernel(radius + i, radius + j) = value;
        }
    }

    double koef = 1 / sum;
    for (int i = 0; i < kernel_size; ++i)
        for (int j = 0; j < kernel_size; ++j)
            kernel(i, j) *= koef;

    return custom(src_image, kernel);
}

/* TODO: separable guassian */

Image gaussian_separable(Image src_image, double sigma, int radius) {
    return src_image;
}

Image median(Image src_image, int radius) {
    return src_image;
}

Image median_linear(Image src_image, int radius) {
    return src_image;
}

Image median_const(Image src_image, int radius) {
    return src_image;
}

static int find_row(Image &img, int row1, int row2) 
{
    int n_rows = img.n_rows;
    int n_cols = img.n_cols;
    
    if (row1 > row2 || row2 >= n_rows)
        return -1;

    int max_sum = 0, max_row = -1;

    for (int i = row1; i <= row2; i++) {
        int sum = 0;

        for (int j = 0; j < n_cols; j++) {
            if (brightness(img(i, j)) >= 1)
                sum++;
        }

        if (sum > max_sum) {
            max_sum = sum;
            max_row = i;
        }
    }

    return max_row;
}

static void suppress_row(Image &img, int row) 
{
    int n_rows = img.n_rows;
    int n_cols = img.n_cols;

    for (int i = max(0, row - kCutMagic); i < min(n_rows, row + kCutMagic); i++) {
        for (int j = 0; j < n_cols; j++) {
            img(i, j) = make_tuple(0, 0, 0);
        }
    }
}

static int find_col(Image &img, int col1, int col2) 
{
    int n_rows = img.n_rows;
    int n_cols = img.n_cols;

    if (col1 > col2 || col2 >= n_cols)
        return -1;

    int max_sum = 0, max_col = -1;

    for (int j = col1; j <= col2; j++) {
        int sum = 0;

        for (int i = 0; i < n_rows; i++) {
            if (brightness(img(i, j)) >= 1)
                sum++;
        }

        if (sum > max_sum) {
            max_sum = sum;
            max_col = j;
        }
    }

    return max_col;
}

static void suppress_col(Image &img, int col) 
{
    int n_rows = img.n_rows;
    int n_cols = img.n_cols;

    for (int j = max(0, col - kCutMagic); j < min(n_cols, col + kCutMagic); j++) {
        for (int i = 0; i < n_rows; i++) {
            img(i, j) = make_tuple(0, 0, 0);
        }
    }
}

void cut_frame(Image src_image, int &x, int &y, int &rows, int &cols) 
{
    int row1, row2;
    int col1, col2;
    int dRow = src_image.n_rows / 10;
    int dCol = src_image.n_cols / 10;

    row1 = find_row(src_image, 0, dRow);
    suppress_row(src_image, row1);
    row1 = max(row1, find_row(src_image, 0, dRow));

    row2 = find_row(src_image, src_image.n_rows - dRow - 1, src_image.n_rows - 1);
    suppress_row(src_image, row2);
    row2 = min(row2, find_row(src_image, src_image.n_rows - dRow - 1, src_image.n_rows - 1));

    col1 = find_col(src_image, 0, dCol);
    suppress_col(src_image, col1);
    col1 = max(col1, find_col(src_image, 0, dCol));

    col2 = find_col(src_image, src_image.n_cols - dCol - 1, src_image.n_cols - 1);
    suppress_col(src_image, col2);
    col2 = min(col2, find_col(src_image, src_image.n_cols - dCol - 1, src_image.n_cols - 1));

    x = row1;
    y = col1;
    rows = row2 - row1 + 1;
    cols = col2 - col1 + 1;

    Image outImage = src_image.deep_copy();

    for (uint i = 0; i < src_image.n_rows; i++) {
        outImage(i, col1) = make_tuple(255, 0, 0);
        outImage(i, col2) = make_tuple(255, 0, 0);
    }
    for (uint j = 0; j < src_image.n_cols; j++) {
        outImage(row1, j) = make_tuple(255, 0, 0);
        outImage(row2, j) = make_tuple(255, 0, 0);            
    }


    save_image(outImage, "/home/alexmihy/mash1718/task1/template/out/rounded.bmp");
}

static void bfs(Matrix<int> &marks, Matrix<int> &was, int x, int y, int n_rows, int n_cols) 
{
    queue<tuple<int, int>> q;

    q.push(make_tuple(x, y));
    was(x, y) = 1;

    while (!q.empty()) {
        int cx, cy;

        tie(cx, cy) = q.front();
        q.pop();

        for (int k = 0; k < 8; k++) {
            int nx, ny;

            nx = cx + kDx[k];
            ny = cy + kDy[k];

            if (nx < 0 || ny < 0 || nx >= n_rows || ny >= n_cols)
                continue;

            if (!marks(nx, ny) || was(nx, ny))
                continue;

            was(nx, ny) = 1;
            q.push(make_tuple(nx, ny));
        }

    }
}

Image canny(Image src_image, int threshold1, int threshold2) {
    int n_rows  = src_image.n_rows;
    int n_cols = src_image.n_cols;

    Image bluredImage;
    Image Kx, Ky;

    bluredImage = gaussian(src_image, 1.4, 2);

    Image grayImg(n_rows, n_cols);

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            uint br;

            br = brightness(bluredImage(i, j));

            grayImg(i, j) = make_tuple(br, br, br);
        }
    }

    Kx = sobel_x(grayImg);
    Ky = sobel_y(grayImg);

    Matrix<int> marks(n_rows, n_cols), was(n_rows, n_cols);
    Image cannyMap(n_rows, n_cols);

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            uint bx, by;
            double phi, b_phi;

            tie(bx, bx, bx) = Kx(i, j);
            tie(by, by, by) = Ky(i, j);

            phi = atan2(by, bx);
            if (phi < 0)
                phi = 2 * M_PI + phi;
            b_phi = phi > M_PI ? phi - M_PI : phi + M_PI;

            int sector1 = phi / (2 * M_PI) * 8;
            int sector2 = b_phi / (2 * M_PI) * 8;

            int G, G1, G2, br;

            G = G_value(Kx, Ky, i, j);
            G1 = G_value(Kx, Ky, i + kDx[sector1], j + kDy[sector1]);
            G2 = G_value(Kx, Ky, i + kDx[sector2], j + kDy[sector2]);
            br = max(min(G, 255), 0);

            marks(i, j) = was(i, j) = 0;
            cannyMap(i, j) = make_tuple(br, br, br);
            if (G < G1 || G < G2)
                cannyMap(i, j) = make_tuple(0, 0, 0);

            if (G < threshold1) {
                cannyMap(i, j) = make_tuple(0, 0, 0);
            } else if (G > threshold2) {
                marks(i, j) = 2;
            } else {
                marks(i, j) = 1;
            }

        }
    }

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            if (marks(i, j) == 2 && !was(i, j)) 
                bfs(marks, was, i, j, n_rows, n_cols);
        }
    }

    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            if (!was(i, j))
                cannyMap(i, j) = make_tuple(0, 0, 0);
        }
    }


    save_image(cannyMap, "/home/alexmihy/mash1718/task1/template/out/map.bmp");

    return cannyMap;
}
