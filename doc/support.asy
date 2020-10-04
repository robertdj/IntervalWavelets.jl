unitsize(0.75cm, 0.4cm);

draw((0,0) -- (16,0));

for(int i = 0; i < 17; ++i) {
    draw((i, -0.5) -- (i, 0.5));
}

label("0", (0, 0.5), N);
label("1", (16, 0.5), N);
for(int i = 1; i < 8; ++i) {
    label(string(i)+"/8", (2*i, 0.5), N);
}

int kstart = -4;
int kend = 2;
for(int n = 1; n <= 8; ++n) {
    label(string(n), (-0.5, -2*n), W);
    label(string(n), (16.5, -2*n), E);
    draw((0, -2*n) -- (16, -2*n));

    kstart += 2;
    kend += 2;
    for(int k = max(kstart, 0); k <= min(kend, 16); ++k) {
        draw((k, -2*n - 0.5) -- (k, -2*n + 0.5));
    }
}

