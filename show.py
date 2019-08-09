import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sg
import simpleaudio as sa


def plot(data, fs, rate, w_ofs, nperseg=1024, at=180):
    plt.subplot(2, 3, 1 + w_ofs)
    plt.plot(data)
    plt.title('Waveform')

    plt.subplot(2, 3, 2 + w_ofs)

    plt.subplots_adjust(wspace=0.8, hspace=0.6)

    f, t, Zxx = sg.stft(data,
                        fs / rate,
                        nperseg=nperseg,
                        window=sg.chebwin(nperseg, at))

    spectrum = 20 * np.log10(2 * np.abs(Zxx) + np.random.rand(Zxx.shape[0], Zxx.shape[1]) * 1e-12)

    mmax = np.max(spectrum) + 6
    plt.pcolormesh(t, f, spectrum, cmap='jet')
    plt.clim(vmin=mmax-at, vmax=mmax)
    plt.colorbar()
    plt.title('STFT Magnitude')
    plt.ylabel('Frequency [Hz]')
    plt.xlabel('Time [sec]')

    plt.subplot(2, 3, 3 + w_ofs)
    plt.plot(spectrum)
    plt.ylim(mmax-at, mmax)
    plt.title('STFT Magnitude')

if __name__ == '__main__':
    plt.clf()

    data = np.loadtxt('data.csv', delimiter=',')

    fs = data[0]
    rate = data[1]
    fs_rs = fs / rate

    sigl = data[2::2]
    sigr = data[3::2]

    plot(sigl, fs, rate, 0, at=140)
    plot(sigr, fs, rate, 3, at=140)

    plt.show()

    audio_array = np.copy(sigl)
    audio_array *= 32767 / max(abs(audio_array))
    audio_array = audio_array.astype(np.int16)

    play_obj = sa.play_buffer(audio_array, 1, 2, int(fs_rs))
    play_obj.wait_done()
