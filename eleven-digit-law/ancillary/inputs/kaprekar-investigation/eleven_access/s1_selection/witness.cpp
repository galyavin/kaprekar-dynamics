#define S1_NO_MAIN
#include "model.cpp"

int main(int argc, char **argv) {
    if (argc!=4) return 2;
    const U row = std::stoull(argv[1]);
    int entrance;
    const Point p = model_input(indexed_ranks(14078212001ULL,row),entrance);
    if (entrance<0) return 1;
    std::cout << "row " << row << " leading-entry " << entrance << " x";
    for (W n : p.n) std::cout << ' ' << numstr(n);
    std::cout << " / " << numstr(p.den) << '\n';
    for (int i=2; i<4; ++i) {
        const Z base = 192000+std::stoll(argv[i]);
        Rng rng{indexed_seed(14078212007ULL,row)};
        const Model a = model(p,base,rng);
        std::cout << "base " << base << " status " << a.status << " clock " << a.time << '\n';
        int stage = 0;
        for (const auto &s : a.passages) {
            std::cout << "stage " << ++stage << " elapsed " << s.time << " length " << s.length
                      << " guard " << s.guard.row << ' ' << s.guard.b << ' ' << s.guard.c
                      << " decrement " << s.decrement << " lane " << s.residue << " / " << s.index
                      << " lost " << s.lost << '\n';
            std::cout << "incoming " << s.incoming << '\n';
            for (const V &v : s.columns) std::cout << "  column" << v << '\n';
            std::cout << "outgoing " << s.outgoing << '\n';
            for (const V &v : s.next_columns) std::cout << "  column" << v << '\n';
            std::cout << "leading section";
            for (W n : s.after.n) std::cout << ' ' << numstr(n);
            std::cout << " / " << numstr(s.after.den) << '\n';
        }
    }
}
