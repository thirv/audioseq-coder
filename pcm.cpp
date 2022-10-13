#include "pcm.hpp"
using namespace std;


pcm_header get_header(unsigned int frame_size, ifstream &f)
{
    pcm_header head;
    head.bytes.resize(44); //wav-related magic number
    f.read(reinterpret_cast<char*>(&head.bytes[0]), head.bytes.size());
    head.nro_ch = head.bytes[23] << 8 | head.bytes[22];
    head.fs = head.bytes[27] << 24 | head.bytes[26] << 16 |
              head.bytes[25] << 8 | head.bytes[24];
    head.nro_bytes = (head.bytes[7] << 24 | head.bytes[6] << 16 |
                      head.bytes[5] << 8 | head.bytes[4]) - 36;
    auto tmp = head.nro_bytes / (2.0 * head.nro_ch * frame_size); // 16 bit only?
    head.nro_fr = (unsigned int) ceil(tmp);
    head.last_size = ((unsigned int) tmp - (head.nro_fr - 1)) * frame_size;
    return head;
}


vector<uint8_t> header_bytes(unsigned int nro_ch, unsigned int fs, unsigned int nro_frames, unsigned int frame_size)
{
    // so far for 16 bit linear quant
    unsigned int nro_bytes = 2 * nro_ch * nro_frames * frame_size;
    vector<uint8_t> bytes = {'R','I','F','F','\0','\0','\0','\0','W','A','V','E','f','m','t',' ',
                             '\x10','\0','\0','\0','\x01','\0','\x01','\0','\0','\0','\0','\0', // size of wave bytes?
                             '\0','\0','\0','\0','\0','\0','\x10','\0','d','a','t','a','\0','\0','\0','\0'};
    unsigned int tmp = nro_bytes + 36;
    bytes[4] = tmp & 0xFF;
    bytes[5] = tmp >> 8 & 0xFF;
    bytes[6] = tmp >> 16 & 0xFF;
    bytes[7] = tmp >> 24 & 0xFF;
    bytes[22] = nro_ch & 0xFF;
    bytes[23] = nro_ch >> 8 & 0xFF;
    bytes[24] = fs & 0xFF;
    bytes[25] = fs >> 8 & 0xFF;
    bytes[26] = fs >> 16 & 0xFF;
    bytes[27] = fs >> 24 & 0xFF;
    tmp = 2 * nro_ch * fs;
    bytes[28] = tmp & 0xFF;
    bytes[29] = tmp >> 8 & 0xFF;
    bytes[30] = tmp >> 16 & 0xFF;
    bytes[31] = tmp >> 24 & 0xFF;
    tmp = 2 * nro_ch;
    bytes[32] = tmp & 0xFF;
    bytes[33] = tmp>> 8 & 0xFF;
    bytes[40] = nro_bytes & 0xFF;
    bytes[41] = nro_bytes >> 8 & 0xFF;
    bytes[42] = nro_bytes >> 16 & 0xFF;
    bytes[43] = nro_bytes >> 24 & 0xFF;
    return bytes;
}


vector<vector<float>> get_pcm(unsigned int nro_ch, unsigned int frame_size, ifstream &f)
{
    vector<uint8_t> bytes(2*nro_ch*frame_size); // for 16 bit wav
    f.read(reinterpret_cast<char*>(&bytes[0]), bytes.size());
    vector<vector<float>> pcm(nro_ch, vector<float> (frame_size));
    for (unsigned int i=0; i<nro_ch; i++) {
        for (unsigned int j=0; j<frame_size; j++) {
            auto ind = 2 * (j * nro_ch + i);
            auto tmp = (int) (bytes[ind+1] << 8 | bytes[ind]);
            if(tmp>32767)
                tmp = tmp - 65536;
            pcm[i][j] = tmp / 32767.0f;
        }
    }
    return pcm;
}


void put_pcm(const vector<vector<float>> &pcm, unsigned int nro_ch, unsigned int frame_size, ofstream &f)
{
    if(frame_size==0)
        return;
    vector<uint8_t> bytes(2*nro_ch*frame_size); // for 16 bit wav;
    
    for (unsigned int i=0; i<nro_ch; i++) {
        for (unsigned int j=0; j<frame_size; j++) {
            auto tmp = (int) round(pcm[i][j]*32767.0f);
            if (tmp<0)
                tmp = tmp + 65536;
            auto ind = 2 * (j * nro_ch + i);
            bytes[ind] = (uint8_t) (tmp & 0xFF);
            bytes[ind+1] = (uint8_t) (tmp >> 8);
        }
    }
    f.write(reinterpret_cast<char*>(&bytes[0]), bytes.size());
}
