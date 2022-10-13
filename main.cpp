#include "pcm.hpp"
#include "coder.hpp"
#include "filebs.hpp"
#include <iostream>
#include <sstream>
using namespace std;

#define FRAME_SIZE 1024
#define FS 48000 // sampling freq


struct sequence{
    vector<settings> s;
    vector<unsigned int> offsets; //bytes from beginning of bs file
    vector<unsigned int> start_frames; //frames from beginning of out file
    vector<unsigned int> end_frames; //frames from beginning of out file
    vector<unsigned int> obj_ids;
    unsigned int num_frames;
    unsigned int num_stems;
    unsigned int fr;
};


sequence init_seq(unsigned int num_objs, vector<uint32_t> nums_frames, vector<uint8_t> obj_rates,
                  vector<vector<float>> start_times, unsigned int frame_size, unsigned int fs)
{ //todo: multichannel objects
    sequence seq;
    seq.num_frames = 0;
    unsigned int offset = 1 + num_objs * 4; // sizeof(uint32_t)
    for (unsigned int obj=0; obj<num_objs; obj++) {
        for (unsigned int objstem=0; objstem<start_times[obj].size(); objstem++) {
            seq.s.push_back(init_coder(frame_size, fs)); // todo: rate input from (unsigned int)obj_rates
            reset_outbuff(1, seq.s.back()); //todo: >1 output channels?
            seq.offsets.push_back(offset);
            seq.start_frames.push_back(round(start_times[obj][objstem]*fs/frame_size)); // rounded to beginning of frame
            seq.end_frames.push_back(seq.start_frames.back()+nums_frames[obj]);
            seq.obj_ids.push_back(obj);
            seq.num_frames = max(seq.num_frames, seq.end_frames.back());
        }
        offset += nums_frames[obj] * seq.s.back().bytes_fr;
    }
    seq.num_stems = (unsigned int) seq.s.size();
    seq.fr = 0;
    return seq;
}


int main(int argv, char *argc[])
{
    //encode todo: output num channels, rates and pannings from input file
	ifstream seq_in("seq.txt",ios::in);
    if (seq_in.fail()) {	
        cerr << "Error: " << strerror(errno) << endl;
        return EXIT_FAILURE;
    }
    string str;
    ofstream f_bsout("bs.bin", ios::out | ios::binary);
    vector<string> file_names;
    vector<vector<float>> start_times;
    vector<uint32_t> nums_frames;
    vector<uint8_t> obj_rates;
    while (getline(seq_in,str)) {
        istringstream in(str);
        string file_name;
        in >> file_name;
        file_names.push_back(file_name);
        vector<float> st;
        float value;
        while(in >> value)
            st.push_back(value);
        start_times.push_back(st);
        ifstream f_in(file_name, ios::in | ios::binary);
        auto head = get_header(FRAME_SIZE, f_in); //testing: auto head2 = header_bytes(1, FS, head.nro_fr, FRAME_SIZE);//579620
        nums_frames.push_back(head.nro_fr-1); // last frame not whole so as tmp hack we discard it
        f_in.close();
        if (head.nro_ch!=1)
            return EXIT_FAILURE; //todo: multichannel objects
    }
    f_bsout.put(static_cast<char>(start_times.size())); // num of audio objects
    f_bsout.write(reinterpret_cast<char*>(&nums_frames[0]), nums_frames.size()*sizeof(nums_frames[0]));
    f_bsout.write(reinterpret_cast<char*>(&obj_rates[0]), obj_rates.size()*sizeof(obj_rates[0]));
    for (unsigned int obj=0; obj<start_times.size(); obj++) {
        ifstream f_in(file_names[obj], ios::in | ios::binary);
        f_in.seekg(44, ios_base::beg); // wav header not needed
        auto s = init_coder(FRAME_SIZE, FS);
        reset_inbuff(1, s);
        obj_rates.push_back((uint8_t)(s.target_rate/1000.0));
        while (s.fr<nums_frames[obj]) {
            auto pcm = get_pcm(1, FRAME_SIZE, f_in);
            auto fr_bs = encoder(pcm,s);
            write_frbs(fr_bs, f_bsout);
            s.fr += 1;
        }
       f_in.close();
    }
    f_bsout.close();
    seq_in.close();
    
    //decode
    ifstream f_bsin("bs.bin", ios::in | ios::binary);
    ofstream f_wavout("out.wav", ios::out | ios::binary);
    char num_objs_char;
    f_bsin.get(num_objs_char);
    auto num_objs = static_cast<unsigned int>(num_objs_char);
    vector<uint32_t> nums_frames_dec(num_objs);
    f_bsin.read(reinterpret_cast<char*>(&nums_frames_dec[0]), nums_frames_dec.size()*sizeof(nums_frames_dec[0]));
    vector<uint8_t> obj_rates_dec(num_objs);
    f_bsin.read(reinterpret_cast<char*>(&obj_rates_dec[0]), obj_rates_dec.size()*sizeof(obj_rates_dec[0]));
    auto seq = init_seq(num_objs, nums_frames_dec, obj_rates_dec, start_times, FRAME_SIZE, FS);
    auto head_bytes  = header_bytes(1, FS, seq.num_frames, FRAME_SIZE);
    f_wavout.write(reinterpret_cast<char*>(&head_bytes[0]), head_bytes.size());
    while (seq.fr<seq.num_frames) {
        vector<vector<float>> pcm_out(1,vector<float>(FRAME_SIZE,0)); //todo: >1 output channels?
        for (unsigned int stem=0; stem<seq.num_stems; stem++) {
            if ((seq.fr>=seq.start_frames[stem]) & (seq.fr<seq.end_frames[stem])) {
                f_bsin.seekg(seq.offsets[stem]+seq.s[stem].fr*seq.s[stem].bytes_fr, ios_base::beg);
                auto fr_bs = read_frbs(f_bsin, seq.s[stem].bytes_fr);
                auto pcm_stem = decoder(fr_bs, seq.s[stem]);
                seq.s[stem].fr += 1;
                transform(pcm_out[0].begin(),pcm_out[0].end(),pcm_stem[0].begin(),pcm_out[0].begin(),plus<float>()); //todo: DRC
            }
        }
        if(seq.fr>1)
            put_pcm(pcm_out, 1, FRAME_SIZE, f_wavout);
        seq.fr += 1;
    }
    f_bsin.close();
    f_wavout.close();

	// todo: stand alone enc dec modes?

	return 0;
}