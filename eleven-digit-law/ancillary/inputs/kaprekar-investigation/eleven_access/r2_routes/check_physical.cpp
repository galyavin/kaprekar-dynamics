#define R2_NO_MAIN
#include "diagnostic.cpp"

int main() {
    Rng rng{13634262519};
    std::vector<Q> path;
    int checked = 0;
    for (int base : {2, 5, 17, 101, 997, 12003, 192015}) {
        for (int row = 0; row < 24; ++row) {
            Word initial{};
            for (int &digit : initial) digit = int(rng.below(base));
            const Physical actual = observe(initial,base,path);
            std::map<Word,int> seen;
            Word w = initial;
            int t = 0;
            for (; !seen.count(w); ++t) {
                seen[w] = t;
                const Q q = word_differences(w);
                w = schoolbook(q,base);
                assert(physical_phase(q,base).next==word_differences(w));
                ++checked;
            }
            assert(actual.clock==seen.at(w));
        }
    }
    std::cout << "168 complete independent schoolbook clocks and " << checked
              << " individual difference transitions agreed\n";
}
