import matplotlib.pyplot as plt
import numpy as np
import math

import os

canal_GSM = 200e3
inceput_banda_GSM = 925e6
sfarsit_banda_GSM = 960e6
banda_GSM = sfarsit_banda_GSM - inceput_banda_GSM
# fft_size = 4096
# nr_esant_sec = 80
# frecventa_start = 924e6
# frecventa_stop = 961.2e6
# fisier = "test80_4096_(924-961.2).dat"
# nr_sec = int(np.ceil(os.path.getsize("test2_80_4096_(924-961.2).dat") / fft_size / nr_esant_sec / 4)) #datele sunt scrise pe 4 byte (sunt float 32)
# pas = fft_size * canal_GSM / (frecventa_stop - frecventa_start)#atunci cand voi parcurge din 200 in 200 kHz

while True:
    #################### Pregatire set date potrivit
    print("\nIntroduceti frecventa de start (multiplu de 200 kHz): ")
    frecventa_start = int(float(input()))
    if frecventa_start % canal_GSM != 0:
        print("Nu ati introdus o frecventa de start potrivita. Trebuie sa fie multiplu de 200kHz")
        continue

    print("Introduceti frecventa de stop: ")
    frecventa_stop = int(float(input()))

    print("Introduceti fft_size: ")
    fft_size = int(float(input()))
    # print("Introduceti cate esantioane pe secunda sunt analizate")
    # nr_esant_sec = int(float(input()))
    pas = fft_size * canal_GSM / (frecventa_stop - frecventa_start)  # atunci cand voi parcurge din 200 in 200 kHz

    prima_decimala = int((pas - math.floor(
        pas)) * 10)  # Aproximez prin lipsa, pentru ca asta ma va obliga sa cresc banda analizata, nu sa o scad si sa imi strice analiza

    if prima_decimala == 0:  # consider ca daca prima zecimala este 0, eroare ce va aparea va fi neglijabila
        pas = round(pas)
        print("Pasul este: " + str(pas))

    else:
        banda_impusa = fft_size * canal_GSM / round(pas) / 1e6
        print("Modificati frecventa de start si cea de stop astfel incat sa aveti banda analizata de: " + str(
            banda_impusa)[:4] + " MHz")
        continue

    print("Introduceti numele fisierului in care sunt datele de analizat")
    fisier = input()

    ############ Citire, dispunere in matrice si calcularea anvelopei cu maxime
    fisier = open(fisier, "rb")  # fisierul este de tip binar
    continut = np.fromfile(fisier, dtype=np.float32)  # datele in gnu sunt salvate in format float32
    nr_linii = int(len(continut) / fft_size)
    nr_coloane = fft_size
    matrice_perioade = np.zeros((nr_linii, nr_coloane))
    anvelopa_maxime = np.zeros(nr_coloane)
    k = 0  # pozitie cursor in continut

    for i in range(nr_linii):  # Dispunerea fisierului pe care il citesc intr-o matrice
        for j in range(nr_coloane):  # Fiecare linie va avea cate o captura fft
            matrice_perioade[i][j] = continut[k]
            k = k + 1

    anvelopa_maxime = [max(i) for i in np.transpose(
        matrice_perioade)]  # Calcularea maximelor pentru a gasi anvelopa care ne va ajuta sa vedem mai usor canalele GSM

    ########### Determinare prag zgomot
    # Stim ca in afara benzii de 925-960 nu se emite nimic, deci e zgomot
    banda_zgomot_start = inceput_banda_GSM - frecventa_start  # banda de zgomot la inceput
    banda_zgomot_stop = frecventa_stop - sfarsit_banda_GSM  # banda de zgomot la final
    banda_frecvente_zgomot_max = np.append(
        anvelopa_maxime[0:int(fft_size * banda_zgomot_start / (frecventa_stop - frecventa_start))],
        anvelopa_maxime[fft_size - int(fft_size * banda_zgomot_stop / (frecventa_stop - frecventa_start)):fft_size])
    nivel_zgomot_max = int(np.mean(banda_frecvente_zgomot_max))

    print("Pragul maxim de zgomot este: " + str(nivel_zgomot_max) + " dB")

    ########### Determinare frecvente centrale
    abatere = 4  # cati dB permitem abaterea fata de o anumita valoare
    esantion_pas = 0  # porneste de la inceputul benzii cu un pas determinat si incearca sa gaseasca doar frecventele centrale ale canalelor GSM
    esantion_pas_imp = int(pas / 2)  # daca pasul merge din valori pare, pas/2 imi va porni de la o valoare impara
    matrice_valori_frecvente_centrale = np.zeros([4, int(fft_size / pas)])

    # pe prima linie imi arata frecventele centrale
    # pe a doua frecventele centrale + 100kHz
    # pe a treia daca sunt GSM sau nu
    # pe a patra va fi un prag ce va folosi la identificarea nivelului mediu de semnal pe canal

    nr_esant_supliment = int(fft_size * 80e3 / (
                frecventa_stop - frecventa_start))  # consideram ca 80kHz sunt suficienti pentru a face media frecv centrale
    temp_par = np.zeros(nr_esant_supliment)
    temp_imp = np.zeros(nr_esant_supliment)

    for i in range(int(fft_size / pas)):  # ajustez valorile frecv centrale(care nu sunt chiar centrale)
        for j in range(nr_esant_supliment):
            temp_par[j] = anvelopa_maxime[esantion_pas + j]
            temp_imp[j] = anvelopa_maxime[esantion_pas_imp + j]

        matrice_valori_frecvente_centrale[0][i] = np.mean(temp_par)
        matrice_valori_frecvente_centrale[1][i] = min(temp_imp)  # cel care este cu 100kHz in fata matrice_valori_frecv_centrale[0]

        matrice_valori_frecvente_centrale[3][i] = int(
            nivel_zgomot_max + (matrice_valori_frecvente_centrale[0][i] - nivel_zgomot_max) / 4)

        esantion_pas += pas
        esantion_pas_imp += pas

    j = 0

    ################ Identificare tip canal
    i = 0
    while i < int(fft_size / pas):
        if matrice_valori_frecvente_centrale[0][i] < max(banda_frecvente_zgomot_max):  # daca e zgomot
            matrice_valori_frecvente_centrale[2][i] = 1
        else:  # daca e semnal
            if ((matrice_valori_frecvente_centrale[0][i] > (matrice_valori_frecvente_centrale[0][i + 1] - abatere)) and
                matrice_valori_frecvente_centrale[0][i] < matrice_valori_frecvente_centrale[0][i + 1]) or (
                    matrice_valori_frecvente_centrale[0][i] < (
                    matrice_valori_frecvente_centrale[0][i + 1] + abatere) and matrice_valori_frecvente_centrale[0][i] >
                    matrice_valori_frecvente_centrale[0][i + 1]):

                # daca 2 frecv centrale sunt relativ apropiate, ori sunt GSM salt in frecv, ori 4G
                j = i
                nr_aparitii = 1

                while ((matrice_valori_frecvente_centrale[0][j] > (
                        matrice_valori_frecvente_centrale[0][j + 1] - abatere)) and
                       matrice_valori_frecvente_centrale[0][j] < matrice_valori_frecvente_centrale[0][j + 1]) or (
                        matrice_valori_frecvente_centrale[0][j] < (
                        matrice_valori_frecvente_centrale[0][j + 1] + abatere) and matrice_valori_frecvente_centrale[0][
                            j] > matrice_valori_frecvente_centrale[0][j + 1]):

                    if matrice_valori_frecvente_centrale[0][j] > matrice_valori_frecvente_centrale[1][j] + abatere:
                        matrice_valori_frecvente_centrale[2][j] = 2  # canal GSM cu salt frecv

                    else:
                        matrice_valori_frecvente_centrale[2][j] = 3  # canal 4G
                    j += 1

                    if j == int(fft_size / pas) - 1:
                        break

                    nr_aparitii += 1  # pentru a vedea cate canale consecutive respecta aceasta conditie de a avea nivel de semnal apropiat; exista posibilitatea de a fi 2 sau 3 cu nivel asemanator si sa fie doar GSM
                # imi iese din bucla cand valoarea de pe frecventa centrala urmamtoare nu este apropiata=>ultima valoare nu o ia in considerare

                if matrice_valori_frecvente_centrale[2][j - 1] == 3:
                    matrice_valori_frecvente_centrale[2][j] = 3

                elif matrice_valori_frecvente_centrale[2][j - 1] == 2:
                    matrice_valori_frecvente_centrale[2][j] = 2

                if nr_aparitii < 5:  # Daca au fost doar 4 canale la rand cu nivel asemanator de semnal, cel mai probabil sunt GSM
                    for k in range(nr_aparitii):
                        matrice_valori_frecvente_centrale[2][j - k] = 7  # modific in urma unde s-a crezut a fi altceva

                else:
                    i = j  # daca nu pun conditia asta, in momentul in care determina x canale consecutive de altceva, momentul in care va ajunge la iteratia x-iteratie<nr_aparitii imi va citi eronat

            else:  # daca nu sunt apropiate, atunci sunt GSM
                matrice_valori_frecvente_centrale[2][i] = 7  # canal GSM
        i += 1

    ############## Identificare frecventa baliza si valoare medie semnal
    valori_medii = np.zeros([3, int(fft_size / pas)])
    # prima linie pt valorile medii
    # a doua pentru a identifica daca e baliza sau nu
    # a treia arata contorul
    contor = 0
    sum = 0
    temp = np.zeros(int(nr_esant_supliment / 2))

    for c in range(int(fft_size / pas)):
        for j in range(nr_linii):
            for k in range(
                    int(nr_esant_supliment / 2)):  # pentru diferentierea tipurilor de canale, media pe mai multe esantioane era eficienta. pentru nivelul de semnal de pe o anumita frecv trebuie sa fim mai rigurosi
                temp[k] = matrice_perioade[j][c * pas + k]

            if matrice_valori_frecvente_centrale[0][c] < max(
                    banda_frecvente_zgomot_max):  # daca e zgomot permanent (Evit cazul in care pe o frecv nu se emite temporar si am zgomot)
                valori_medii[0][c] = np.mean(temp)
                break  # daca stiu sigur ca e zgomot, imi e de ajuns sa verific o singura linie, nu pe toate

            else:

                if max(temp) > matrice_valori_frecvente_centrale[3][c]:  # iau in considerare doar valorile peste prag
                    contor += 1
                    sum += np.mean(temp)
                    valori_medii[0][c] = sum / contor
                    valori_medii[2][c] = contor

        if contor > int(nr_linii * 0.85) and matrice_valori_frecvente_centrale[0][c] > nivel_zgomot_max:
            valori_medii[1][c] = 1  # am baliza

        sum = 0
        contor = 0

    ################ Afisari

    for i in range(int(fft_size / pas)):
        if matrice_valori_frecvente_centrale[2][i] == 2:
            print("La frecventa " + str((frecventa_start + i * 200e3) / 1e6)[:5] + " MHz, avem GSM Freq Hopping si nivel mediu de semnal: " + str(
                (valori_medii[0][i]))[:5] + " dB")

        elif matrice_valori_frecvente_centrale[2][i] == 3:
            print("La frecventa " + str((frecventa_start + i * 200e3) / 1e6)[:5] + " MHz, avem 4G si nivel mediu de semnal: " + str((valori_medii[0][i]))[:5] + " dB")

        elif matrice_valori_frecvente_centrale[2][i] == 7:
            if valori_medii[1][i] == 1:
                print("La frecventa " + str((frecventa_start + i * 200e3) / 1e6)[:5] + " MHz, avem GSM - BALIZA si nivel mediu de semnal: " + str((valori_medii[0][i]))[:5] + " dB")
            else:
                print("La frecventa " + str((frecventa_start + i * 200e3) / 1e6)[:5] + " MHz, avem GSM - TRAFIC si nivel mediu de semnal: " + str((valori_medii[0][i]))[:5] + " dB")
        elif matrice_valori_frecvente_centrale[2][i] == 1:
            print("La frecventa " + str((frecventa_start + i * 200e3) / 1e6)[:5] + " MHz, avem zgomot si nivel mediu de semnal: " + str((valori_medii[0][i]))[:5] + " dB")
    frecvente_centrale = np.linspace(frecventa_start, frecventa_stop, int(fft_size / pas))
    axa_frecvente = np.linspace(frecventa_start, frecventa_stop, fft_size)

    plt.figure()
    plt.plot(axa_frecvente, anvelopa_maxime, 'r')
    plt.plot(axa_frecvente, continut[0:fft_size], 'b')

    plt.figure()
    plt.stem(frecvente_centrale, matrice_valori_frecvente_centrale[0])
    plt.title("Valori frecvente centrale")

    plt.figure()
    plt.stem(frecvente_centrale, matrice_valori_frecvente_centrale[1], 'black')
    plt.title("Valori frecvente centrale + 100 kHz")
    plt.show()