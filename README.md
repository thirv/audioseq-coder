# audioseq-coder
Compact audio coder for audio sequences. Motivation is to use pure modern C++ (currently C++11) limited to standard library. Somewhat functional programming style, which suits audio processing in my oppinion better than OO.

Unlike typical audio codecs that operate on audio frames, this is tailored toward compressing separate audio files and repeating them according to the pattern playing defined in seq.txt. This is suitable for music that has many audio elements/segments/objects (short or long) that repeat, such as electronic music. The principles used in the compression are typical ones used in audio compression; MDCT transform, quantization, and entropy coding to a maximally compact bitsream. Some preliminary ideas for exploiting correlations between the sequences/objects via PCA (https://en.wikipedia.org/wiki/Principal_component_analysis), and compressing only the largest eigenvectors.

- main.cpp is the entry point for a small running example.
- coder.cpp includes the encoder and decoder to generate and extract audio bitsream. 
- bs.cpp has the entropy coding
- filebs.cpp includes bitsream construction via manipulating bytes and file writing
- pcm.hpp implements audio file reading
- mdct.cpp implementation of https://en.wikipedia.org/wiki/Modified_discrete_cosine_transform
- pvq.cpp implementation of https://en.wikipedia.org/wiki/Pyramid_vector_quantization


TODO: Code is compact with minimal comments, better documentaion. Feel free to import to your favorite IDE etc.
