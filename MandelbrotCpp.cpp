//Maurycy Muzyka 285730
//WdPRiR 2023Z

//Mandelbrot w C++

//https://code.visualstudio.com/docs/cpp/config-mingw

#include <iostream>
#include <complex> //https://en.cppreference.com/w/cpp/numeric/complex
#include <cmath>
#include <chrono>
#include <thread>
#include <vector>
#include <future>
#include <mutex>
#include <execution>
#include <exception>
using namespace std;

double liczZbieznosc(std::complex<double> z0, int maxx) {
    std::complex<double> z = z0;
    for (int t = 0; t < maxx; t++) {
        double fAbs = abs(z);
        if (fAbs > 2.0)
        {
             // based on the final value, add a fractional amount based on
             // how much it escaped by (fAbs will be in the range of 2 to around 4):
             return (float) (t + (2.0 - (log(fAbs) / log(2.0))));
        }
        z = z * z + z0;
    }
    return maxx;
}
	
double** liczBlok(double x_start3, int x_n3, int y_n3, double y_end3, double dx3, double dy3) {
	double** wyn = new double*[x_n3];
	for(int i = 0; i < x_n3; i++) {
        wyn[i] = new double[y_n3];
		for(int j = 0; j < y_n3; j++) {
			std::complex<double> z(x_start3 + dx3*(i + 0.5), y_end3 - dy3*(j + 0.5));
			wyn[i][j] = liczZbieznosc(z, 100);
			//wyn[i][j] = real(z); //debug czy dobrze wycinkuje obszar
		}
	}
	return wyn;
}

int main() { //throws InterruptedException, ExecutionException {
	cout << "Start" << endl;

	double x_start = -2.;
	double x_end = -x_start;
	double y_start = x_start;
	double y_end = x_end;

	int analiza = 1; //1 - rozne rozmiary bloku; 2 - rozne px
	
	int px[] = {1000}; 
	int rozmBloku[] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000}; //rozmiary blokow
    //int px[] = {10};
    //int rozmBloku[] = {1, 2, 5, 10};
	if(analiza == 2) {
		int px[] = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000}; //rozdzielczosci
        int rozmBloku[] = {1};
	}
	int nBlokow[max(end(rozmBloku)-begin(rozmBloku), end(px)-begin(px))]; //liczba blokow
	for(int i = 0; i < end(rozmBloku)-begin(rozmBloku); i++) {
		if(analiza == 1)
			nBlokow[i] = (int)(px[0] / rozmBloku[i]);
		else
			nBlokow[i] = (int)(px[i] / rozmBloku[0]);
	}	
		
	int zapis = 0; //1 - zapisz csv; UWAGA! Wczesniej ustaw 1-elementowy rozmiar bloku i 1-elementowa rozdzielczosc px
	int printuj = 1; //1 - jesli zapisujesz to wyswietl dodatkowo "obraz" w terminalu
	
	int n_usr = 10; 
	for(int k = 0; k < end(px)-begin(px); k++) { //rozna rozdzielczosc obrazka
		double t_sr = 0;
		int x_n = px[k];
		int y_n = x_n;
		double dx = (x_end-x_start)/x_n;
		double dy = (y_end-y_start)/y_n;
		double** wyniki = new double*[x_n];
        for(int i = 0; i < x_n; i++)
            wyniki[i] = new double[y_n];
			
		for(int m = 0; m < end(rozmBloku)-begin(rozmBloku); m++) { //rozne rozmiary blokow
            //cout << "m: " << m << endl;
			int x_n2 = (int)(x_n / nBlokow[m]);
			int y_n2 = y_n; 
			
			for(int l = 0; l < n_usr; l++) { //usrednianie realizacji calego procesu

                std::vector<std::future<double**>> results;
                std::vector<std::thread> threads;
                //LinkedList<Future<double[][]>> results = new LinkedList<>();
				//ExecutorService ex = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

                auto start = std::chrono::high_resolution_clock::now();
                //long start = System.nanoTime();

                for(int o = 0; o < nBlokow[m]; o++) { //kolejne bloki jako joby
                    int p = o;
                    double x_start2 = (double)(x_start + p * dx * rozmBloku[m]);
					double y_end2 = y_end;
                    results.emplace_back(std::async(std::launch::async, liczBlok, x_start2, x_n2, y_n2, y_end2, dx, dy));
                    //results.add(ex.submit(() -> liczBlok(x_start2, x_n2, y_n2, y_end2, dx, dy))); //liczBlok()
                }

                for(int o = 0; o < nBlokow[m]; o++) {
                    int p = o;
                    int i2 = (int)(p * rozmBloku[m]); //przesuniecie x-owe do tablicy wynikow[][]
                    try {
                        auto block = results[p].get();
                        for(size_t j = 0; j < y_n2; j++) {
                            for(size_t i = 0; i < x_n2; i++) {
                                wyniki[i + i2][j] = block[i][j];
                            }
                        }
                    } catch (const std::exception& e) {}
                    /*try {
						for(int j = 0; j < y_n2; j++) {
							for(int i = 0; i < x_n2; i++) {
								wyniki[i + i2][j] = (double)(results.get(p).get()[i][j]);  //czeka aż wynik sie pojawi
							}
						}
					} catch (InterruptedException e) {} catch (ExecutionException e) {}*/
                }
				//ex.shutdown();
				//ex.awaitTermination(1, TimeUnit.DAYS);
					            
                /*
                for(int i = 0; i < x_n; i++) {
					for(int j = 0; j < y_n; j++) {
						std::complex<double> z(x_start + dx*(i + 0.5), y_end - dy*(j + 0.5));
						wyniki[i][j] = liczZbieznosc(z, 100);
					}
				}
                */

                auto stop = std::chrono::high_resolution_clock::now();
                t_sr += std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start).count()/1e9;
                /*long stop = System.nanoTime();
				t_sr += (double)(stop - start)/1e9;*/
			}
			if(analiza == 1) {
				cout << "Rozmiar bloku: " << rozmBloku[m] << endl;
				cout << "Time: " << t_sr/n_usr << endl;
			}
			t_sr = 0;
			
			if(zapis == 1) {
                if(printuj == 1) {    
		            for(int j = 0; j < y_n; j++) {
						for(int i = 0; i < x_n; i++) {
                            cout << wyniki[i][j] << " ";
						}
                    cout << endl;
				    }
		        }
			}
		}
		if(analiza == 2) {
			cout << "Pixeli: " << x_n*y_n << endl;
            cout << "Time: " << t_sr/n_usr << endl;
		}
	}
    cout << "Koniec" << endl;
    return 0;
}