#include <Siv3D.hpp>

#include <numbers>
#include <algorithm>
#include <type_traits>
#include <utility>
#include <ppl.h>

std::vector<double> lpc(const Array<WaveSampleS16>& input, std::size_t posSample);
std::vector<double> freqz(const std::vector<double>& b, const std::vector<double>& a);

void Main() {
    // マイクをセットアップ
    Microphone mic(unspecified);

    if(!mic) {
        // マイクを利用できない場合、終了
        throw Error(U"Microphone not available");
    }

    // 録音をスタート
    mic.start();

    FFTResult fft;

    while(System::Update()) {
        // FFT の結果を取得
        mic.fft(fft);
		const Array<WaveSampleS16>& buf = mic.getBuffer();
		const std::size_t posSample = mic.posSample();
		// std::vector<double> t = lpc(buf, 64, posSample, fft.samplingRate / buf.size());
		std::vector<double> t = lpc(buf, posSample);
		double mul = *std::max_element(fft.buffer.begin(), fft.buffer.end()) / *std::max_element(t.begin(), t.end());
		for (auto& te : t) te *= mul;

        // 結果を可視化
		int wid = 8192 / (8192);

		std::vector<std::pair<int, int>> formantPosList;
		for (auto i : step(4096)) {
			if (0 < i && i + 1 < t.size() && t[i - 1] < t[i] && t[i] >= t[i + 1]) {
				int size2 = Pow(t[i], 0.6f) * 1200 * 2 * 2;
				formantPosList.push_back({ i, size2 });
			}
		}

        for(auto i : step(4096)) {
			const double size = Pow(fft.buffer[i], 0.6f) * 1200 * 2;
			double size2 = 0;

			if (t.size() > i)
				size2 = Pow(t[i], 0.6f) * 1200 * 2 * 2;
            // const double size = (Log(fft.buffer[i])/Log(10) + 5) * 150;

			if (0 < i && i+1 < t.size() && t[i-1] < t[i] && t[i] >= t[i+1])
				RectF(Arg::bottomLeft(i * wid, 600), wid, 600).draw((ColorF(1)));

			auto col = HSV(240 - i);
			col.a = 0.4;
            RectF(Arg::bottomLeft(i, 600), 1, size).draw(col);
			if(t.size() > i)
				RectF(Arg::bottomLeft(i * wid, 600-size2), wid, 2).draw((ColorF(0.5) + HSV(240 - i)));
        }

		for(auto i : step(formantPosList.size()-1))
			Line(
				formantPosList[i  ].first, 600 - formantPosList[i  ].second,
				formantPosList[i+1].first, 600 - formantPosList[i+1].second
			).draw(LineStyle::RoundCap, 2, Palette::Orange);

        // 周波数表示
        Rect(Cursor::Pos().x, 0, 1, Scene::Height()).draw();
        ClearPrint();
		// Print << U"{} Hz \nfft: {} / buf: {}"_fmt(Cursor::Pos().x * fft.resolution, fft.buffer.size(), buf.size());
		auto f_0 = formantPosList.empty() ? 0 : formantPosList[0].first;
		Print << U"f_0   : {} Hz\ncursor: {} Hz"_fmt(f_0 * fft.resolution, Cursor::Pos().x * fft.resolution);
    }
}

// TODO: バッファーをもう少しまともな場所に置く
constexpr int sampleSize = 8192;
constexpr int order = 128;
constexpr int lags_num = order + 1;
std::vector<double> a(order + 1), e(order + 1);
std::vector<double> R(order + 1);  // 自己相関関数
std::vector<double> U(order + 3);
std::vector<double> V(order + 3);
std::vector<double> H(sampleSize + 1);

std::vector<double> lpc(const Array<WaveSampleS16>& input, std::size_t posSample) {
	assert(sampleSize < input.size());
	// バッファーから index を計算する関数 ()
	auto inputAt =
		[=](std::size_t idx) {
			assert(idx <= sampleSize);
			const auto pos =
				posSample >= idx
				? posSample - idx
				: input.size() + posSample - idx;
			const auto& sample = input[pos];
			return
				(static_cast<std::int_fast64_t>(sample.left) +
				 static_cast<std::int_fast64_t>(sample.right)) / 2;
		};
	// 自己相関関数 の計算
	for (int dt : step(order + 1)) {
		R[dt] = 0;
		for (int n : step(sampleSize - dt))
			R[dt] += inputAt(n) * inputAt(n + dt) / 32768.0;
	}
	

	// Levinson Durbin algorithm
	// https://tmytokai.github.io/open-ed/activity/d-filter/text04/page04.html とかが参考になる？
	// 講義資料： https://ahcweb01.naist.jp/lecture/2016/sp/material/sp-v2.pdf

	for (auto& t : a) t = 0;
	for (auto& t : e) t = 0;

	a[0] = 1;
	e[0] = 1;
	a[1] = -R[1] / R[0];
	e[1] =  R[0] + R[1] * a[1];

	for(int k : step(1, order - 1)) {  // [ 1, 1+(order-1) )
		double lambda{};
		for(int j : step(k + 1))
			lambda -= R[k + 1 - j] * a[j];
		lambda /= e[k];

		U[0] = 1.0; V[0] = 0.0;
		for(int i : step(k + 1)) {
			U[i] = a[i];
			V[k + 1 - i] = a[i];
		}
		U[k + 1] = 0.0;
		V[k + 1] = 1.0;

		for(int i : step(k + 2))
			a[i] = U[i] + lambda * V[i];

		e[k + 1] = e[k] * (1.0 - lambda * lambda);
	}

	return freqz(e, a);
}

std::vector<double> freqz(const std::vector<double>& e, const std::vector<double>& a) {
	std::complex<double> zn{};

	for(int n : step(sampleSize + 1)) {
	// Concurrency::parallel_for(0, sampleSize + 1, 1, [&](int n) {
		auto z = std::exp(std::complex<double>(0.0, -2.0 * std::numbers::pi * n / sampleSize));
		std::complex<double> numerator{}, denominator{};

		zn = 1;
		for (auto itr = e.rbegin(); itr != e.rend(); ++itr) {
			numerator += *itr * zn;
			zn *= z;
		}

		zn = 1;
		for (auto itr = a.rbegin(); itr != a.rend(); ++itr) {
			denominator += *itr * zn;
			zn *= z;
		}

		H[n] = std::abs(numerator / denominator);
	};

    return H;
}
