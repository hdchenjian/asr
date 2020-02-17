from python_speech_features import mfcc
from python_speech_features import delta
from python_speech_features import logfbank
import scipy.io.wavfile as wav

from scikits.talkbox.features import mfcc as mfcc_new

(rate, sig) = wav.read("asr/audio/train/00_01_giuseppe.s.wav")
print rate, sig[0:10]
mfcc_feat = mfcc(sig, rate, nfft=512*2, winlen=0.023)
print mfcc_feat.shape
print mfcc_feat[0]

#features = mfcc_new(sound, nwin=, fs=rate, nceps=13)
exit()

d_mfcc_feat = delta(mfcc_feat, 2)
print d_mfcc_feat.shape


fbank_feat = logfbank(sig, rate)

print fbank_feat.shape

#print(fbank_feat[1:3,:])
