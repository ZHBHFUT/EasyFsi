#include <fstream>
#include <string>
#include <iomanip>

#include "ModalResults.hpp"
#include "Logger.hpp"

#ifdef _MSC_VER
#pragma warning(disable:4996)  // fopen
#endif

namespace EasyLib {
    void ModalResults::clear()
    {
        ngrid = nmode = 0;
        ids.clear();
        coords.clear();
        freq.clear();
        mass.clear();
        phi.clear();
    }
    void ModalResults::load(const char* file)
    {
        std::ifstream ifs(file);
        if (!ifs)error("failed open file: %s", file);

        info("reading modal results file: %s\n", file);

        auto skip_comment = [](std::ifstream& ifs) {
            std::string s;
            while (!ifs.eof()) {
                auto pos = ifs.tellg();
                std::getline(ifs, s);
                if (s.empty() || s.front() == '$')continue;
                ifs.seekg(pos);
                break;
            }
        };

        // NGRID NMODE
        skip_comment(ifs);
        ifs >> ngrid >> nmode;
        info(
            "    NGRID = %d\n"
            "    NMODE = %d\n",
            ngrid, nmode
        );

        // allocate
        ids.resize(ngrid, 0);
        coords.resize(ngrid);
        freq.resize(nmode, 0.0);
        mass.resize(nmode, 1.0);
        phi.resize(nmode, ngrid, 6);

        // GRID X Y Z
        skip_comment(ifs);
        for (int i = 0; i < ngrid; ++i)
            ifs >> ids[i] >> coords[i].x >> coords[i].y >> coords[i].z;

        // modal data
        for (int imod = 0; imod < nmode; ++imod) {
            // FREQ  MASS
            skip_comment(ifs);
            ifs >> freq[imod] >> mass[imod];
            // shape
            skip_comment(ifs);
            int id;
            for (int i = 0; i < ngrid; ++i) {
                ifs >> id
                    >> phi(imod, i, 0) >> phi(imod, i, 1) >> phi(imod, i, 2)
                    >> phi(imod, i, 3) >> phi(imod, i, 4) >> phi(imod, i, 5);
                if (id != ids[i])error("grid id not agree!");
            }
        }
        ifs.close();
        info("!!!OK!!!\n");
    }
    void ModalResults::save(const char* file)const
    {
        FILE* fp = fopen(file, "wt");
        if (fp == nullptr) {
            error("failed open file: %s", file);
            return;
        }

        info("save modal results to file: %s\n", file);

        // NGRID NMODE
        //      1234567890123456
        fprintf(fp, "$%7s %8s\n", "NGRID", "NMODE");
        fprintf(fp, "%8d %8d\n", ngrid, nmode);

        // GRID X Y Z
        fprintf(fp, "$%7s %24s %24s %24s\n", "GRID", "X", "Y", "Z");
        for (int i = 0; i < ngrid; ++i) {
            fprintf(fp, "%8d %24.16E %24.16E %24.16E\n",
                ids[i], coords[i].x, coords[i].y, coords[i].z
            );
        }

        // modal
        for (int imod = 0; imod < nmode; ++imod) {
            // FREQ MASS
            fprintf(fp, "$ MODE = %d\n", imod + 1);
            fprintf(fp, "$%23s %24s\n", "FREQ(rad/s)", "MASS");
            fprintf(fp, "%24.16E %24.16E\n", freq[imod], mass[imod]);
            // GRID UX UY UZ RX RY RZ
            fprintf(fp, "$%7s %24s %24s %24s %24s %24s %24s\n",
                "GRID", "UX", "UY", "UZ",
                "ROTX", "ROTY", "ROTZ");
            for (int i = 0; i < ngrid; ++i) {
                fprintf(fp, "%8d %24.16E %24.16E %24.16E %24.16E %24.16E %24.16E\n",
                    ids[i],
                    phi(imod, i, 0), phi(imod, i, 1), phi(imod, i, 2),
                    phi(imod, i, 3), phi(imod, i, 4), phi(imod, i, 5)
                );
            }
        }
        fclose(fp);
        info("!!!OK!!!\n");
    }
}
