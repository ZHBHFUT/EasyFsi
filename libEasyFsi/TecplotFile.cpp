#include <limits>
#include <fstream>
#include <string>
#include <sstream>

#include "Logger.hpp"
#include "TecplotFile.hpp"

namespace EasyLib {

    //! @brief skip white space including end-of-line
    void skip_whitespace_and_comment(std::istream& is, size_t& nline)
    {
        // skip first comment line
        if (!is.eof() && is.peek() == '#') {
            is.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // '\n' will be extracted from stream.
            ++nline;
        }

        // skip whitespace
        while (!is.eof() && std::isspace(is.peek())) {
            // read one character
            if (is.get() == '\n') {
                ++nline;

                // skip comment line
                if (!is.eof() && is.peek() == '#') {
                    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // '\n' will be extracted from stream.
                    ++nline;
                }
            }
        }
    }

    //! @brief Reading a key
    bool read_key(std::istream& is, std::string& key, size_t& nline)
    {
        skip_whitespace_and_comment(is, nline);

        key.clear();
        
        // key must be starting with alpha character
        if (std::isalpha(is.peek())) {
            // read token
            while (!is.eof() && std::isalnum(is.peek()) && is.peek() != '=') {
                key.push_back((char)is.get());
            }
            if (key.empty()) {
                info("***ERROR*** key is empty");
                return false;
            }
            skip_whitespace_and_comment(is, nline);

            // skip '=' for next reading
            if (!is.eof() && is.peek() == '=')is.get();
        }

        // there is no key if key is empty.
        return true;
    }

    //! @brief Reading string value in format: "???"
    bool read_string_value(std::istream& is, std::string& value, size_t& nline)
    {
        skip_whitespace_and_comment(is, nline);
        
        // must be starting with '"'
        if (is.peek() == '\"') {
            // skip start \"
            is.get();

            // read until match closure '"'
            std::getline(is, value, '\"');

            // check: value should not contain end-of-line.
            if (!value.empty() && value.find_first_of('\n') != std::string::npos) {
                info("***ERROR*** string value cross multi-line!");
                return false;
            }

            // continue reading if last is \"
            while (!is.eof() && !value.empty() && value.back() == '\\') {
                // modify last character as '"'
                value.back() = '\"';

                // skip \"
                is.get();

                // read content between "???"
                std::string s;
                std::getline(is, s, '\"');

                // check: value should not contain end-of-line.
                if (!s.empty() && s.find_first_of('\n') != std::string::npos) {
                    info("***ERROR*** string value cross multi-line!");
                    return false;
                }

                // append
                value.append(s);
            }

            skip_whitespace_and_comment(is, nline);

            // skip "," if it exists
            if (!is.eof() && is.peek() == ',')is.get();

            return true;
        }
        else {
            value.clear();
            return false;
        }
    }

    bool read_string_values(std::istream& is, std::vector<std::string>& values, size_t& nline)
    {
        // "v1", "v2", ...
        // "v1"  "v2"  ...
        
        // read first
        if (!read_string_value(is, values.emplace_back(std::string{}), nline))return false;

        // read more
        for (;;) {
            skip_whitespace_and_comment(is, nline);
            if (is.peek() == '\"') {
                if (!read_string_value(is, values.emplace_back(std::string{}), nline))return false;
            }
            else
                break;
        }
        return true;
    }

    bool read_value(std::istream& is, std::string& value, size_t& nline)
    {
        skip_whitespace_and_comment(is, nline);

        // read until match \"
        if (!is.eof() && is.peek() == '\"') {
            return read_string_value(is, value, nline);
        }
        // read until match space or comma
        else {
            value.clear();
            while (!is.eof() && !std::isspace(is.peek()) && is.peek() != ',') {
                value.push_back((char)is.get());
            }
        }

        skip_whitespace_and_comment(is, nline);

        // skip comma
        if (!is.eof() && is.peek() == ',')is.get();

        return true;
    }

    bool read_int_value(std::istream& is, int& value, size_t& nline)
    {
        skip_whitespace_and_comment(is, nline);
        if (is.eof())return false;

        is >> value;
        if (!is.good())return false;

        skip_whitespace_and_comment(is, nline);

        // skip comma
        if (!is.eof() && is.peek() == ',')is.get();
        return true;
    }

    bool read_dbl_value(std::istream& is, double& value, size_t& nline)
    {
        skip_whitespace_and_comment(is, nline);
        if (is.eof())return false;

        is >> value;
        if (!is.good())return false;

        skip_whitespace_and_comment(is, nline);

        // skip comma
        if (!is.eof() && is.peek() == ',')is.get();
        return true;
    }

    bool read_var_loc(std::istream& is, std::vector<TecVarLocation>& var_loc, size_t& nline)
    {
        // ([]=xx,[]=xxx,...)

        skip_whitespace_and_comment(is, nline);
        if (is.eof())return false;

        if (is.get() != '(') {
            info("\n***ERROR*** syntax error for VARLOCTION\n");
            return false;
        }

        // [x-x,xx,xx]=xx

        while (!is.eof()) {
            if (is.get() != '[') {
                info("\n***ERROR*** syntax error for VARLOCTION\n");
                return false;
            }

            // reading content between []
            std::string s;
            std::getline(is, s, ']');
            
            // parse location type: [...]=location

            char c = '\0';
            is >> c;
            if (c != '=') {
                info("\n***ERROR*** syntax error for VARLOCTION, '=' not found\n");
                return false;
            }
            std::string value;
            if (!read_value(is, value, nline))return false; // comma is skipped
            TecVarLocation loc = NODAL;
            if      (::_stricmp(value.c_str(), "NODAL"       ))loc = NODAL;
            else if (::_stricmp(value.c_str(), "CELLCENTERED"))loc = CELLCENTERED;
            else {
                info("\n***ERROR*** invalid VARLOCTION type: %s\n", value.c_str());
                return false;
            }

            // parse id: x-y,xx,...

            std::istringstream iss(s);
            while (!iss.eof()) {
                int id0 = 0, id1 = 0;
                iss >> id0;
                
                size_t j = 0;
                skip_whitespace_and_comment(iss, j);

                id1 = id0;
                if (!iss.eof()) {
                    char ch = '\0';
                    iss >> ch;
                    if (ch == '-') {
                        iss >> id1;
                        skip_whitespace_and_comment(iss, j);
                        if (!iss.eof())iss >> ch; // c = ','
                    }

                    if (ch != ',') {
                        info("\n***ERROR*** syntax error for VARLOCTION\n");
                        return false;
                    }
                }

                // update location
                for (int i = id0 - 1; i < id1; ++i)
                    var_loc.at(i) = loc;
            }

            if (c == ')')return true; // end of block
        }
        return false;
    }

    bool read_var_share_list(std::istream& is, std::vector<int>& var_share_list, int zlast, size_t& nline)
    {
        // ([xx-y,xxx,...]=z1,[]...=z2,...)

        skip_whitespace_and_comment(is, nline);
        if (is.eof())return false;
        if (is.get() != '(') {
            info("\n***ERROR*** syntax error for VARSHARELIST\n");
            return false;
        }

        // [x-x,xx,xx]=xx

        while (!is.eof()) {
            if (is.get() != '[') {
                info("\n***ERROR*** syntax error for VARSHARELIST\n");
                return false;
            }

            // reading content between []
            std::string s;
            std::getline(is, s, ']');

            //skip_whitespace_and_comment(is, nline);

            // parse zone id: [...]=zone

            int z0 = zlast;

            int c = '\0';
            is >> c;
            if (c == '=')is >> z0 >> c;

            // parse var id: [xx-y,xxx,...]
            std::istringstream iss(s);
            while (!iss.eof()) {
                int id0 = 0, id1 = 0;
                iss >> id0;

                size_t j = 0;
                skip_whitespace_and_comment(iss, j);

                id1 = id0;
                if (!iss.eof()) {
                    char ch = '\0';
                    iss >> ch;
                    if (ch == '-') {
                        iss >> id1;
                        skip_whitespace_and_comment(iss, j);
                        if (!iss.eof())iss >> ch; // ch = ','
                    }

                    if (ch != ',') {
                        info("\n***ERROR*** syntax error for VARLOCTION\n");
                        return false;
                    }
                }

                // update
                for (int i = id0 - 1; i < id1; ++i)
                    var_share_list.at(i) = zlast;
            }

            if      (c == ')')return true;
            else if (c != ',') {
                info("\n***ERROR*** syntax error for VARSHARELIST\n");
                return false;
            }
        }
        return false;
    }

    void TecplotFile::read_ascii(const char* file)
    {
        std::ifstream ifs(file);
        if (!ifs.is_open())error("failed opening tecplot file: %s", file);

        info("\nreading ascii tecplot file: %s\n", file);

        size_t iline = 0;
        std::string key, value;

        // read header
        while (!ifs.eof() && read_key(ifs, key, iline)) {
            
            // not a key
            if (key.empty()) {
                error("failed reading key data"); return;
            }

            // TITLE
            if      (::_stricmp(key.c_str(), "TITLE") == 0) {
                if (!read_value(ifs, title_, iline)) { error("failed reading TITLE value"); return; }
                info("  TITLE = \"%s\"\n", title_.c_str());
            }
            // FILETYPE
            else if (::_stricmp(key.c_str(), "FILETYPE") == 0) {
                if (!read_value(ifs, value, iline)) { error("failed reading FILETYPE value"); return; }
                if      (::_stricmp(value.c_str(), "FULL") == 0)type_ = FULL;
                else if (::_stricmp(value.c_str(), "GRID") == 0)type_ = GRID;
                else if (::_stricmp(value.c_str(), "SOLUTION") == 0)type_ = SOLUTION;
                else {
                    error("invalid file type: %s", value.c_str());
                    return;
                }
                info("  TITLE = %s\n", value.c_str());
            }
            // VARIABLES
            else if (::_stricmp(key.c_str(), "VARIABLES") == 0) {
                if (!read_string_values(ifs, variables_, iline)) { error("failed reading VARIABLES value"); return; }
                info("  VARIABLES =");
                for (auto& s : variables_)info(" \"%s\"", s.c_str());
                info("\n");
            }
            // ZONE
            else if (::_stricmp(key.c_str(), "ZONE") == 0) {
                break;
            }
        }

        // read zones

        zones_.clear();
        for (;;) {
            int zlast = 0;
            auto& z = zones_.emplace_back(TecplotZone{});
            info("  reading zone #%d\n", (int)zones_.size());

            z.var_loc.resize(variables_.size(), NODAL);
            z.var_share_list.resize(variables_.size(), zlast);

            // zone header
            while (!ifs.eof() && read_key(ifs, key, iline)) {
                // not a key: data start
                if (key.empty())break;

                // title
                if (::_stricmp(key.c_str(), "T") == 0) {
                    if (!read_string_value(ifs, z.title, iline)) { error("failed reading T value"); return; }
                    info("    TITLE = %s\n", z.title.c_str());
                }
                //ZONETYPE
                else if (::_stricmp(key.c_str(), "ZONETYPE") == 0) {
                    if (!read_string_value(ifs, value, iline)) { error("failed reading ZONETYPE value"); return; }
                    if (::_stricmp(value.c_str(), "ORDERED") == 0)z.type = ORDERED;
                    else if (::_stricmp(value.c_str(), "FELINESEG") == 0)z.type = FELINESEG;
                    else if (::_stricmp(value.c_str(), "FETRIANGLE") == 0)z.type = FETRIANGLE;
                    else if (::_stricmp(value.c_str(), "FEQUADRILATERAL") == 0)z.type = FEQUADRILATERAL;
                    else if (::_stricmp(value.c_str(), "FETETRAHEDRON") == 0)z.type = FETETRAHEDRON;
                    else if (::_stricmp(value.c_str(), "FEBRICK") == 0)z.type = FEBRICK;
                    else if (::_stricmp(value.c_str(), "FEPOLYGON") == 0)z.type = FEPOLYGON;
                    else if (::_stricmp(value.c_str(), "FEPOLYHEDRAL") == 0)z.type = FEPOLYHEDRAL;
                    else { error("invalid zone type: %s", value.c_str()); return; }
                    info("    TITLE = %s\n", value.c_str());
                }
                // F
                else if (::_stricmp(key.c_str(), "F") == 0) {
                    if (!read_value(ifs, value, iline)) { error("failed reading F value"); return; }
                    z.type = FEQUADRILATERAL;
                    z.is_block = false;
                    info("    F = %s\n", value.c_str());
                }
                else if (::_stricmp(key.c_str(), "I") == 0) {
                    if (!read_int_value(ifs, z.i, iline)) { error("failed reading I value"); return; }
                    info("    I = %d\n", z.i);
                }
                else if (::_stricmp(key.c_str(), "J") == 0) {
                    if (!read_int_value(ifs, z.j, iline)) { error("failed reading J value"); return; }
                    info("    J = %d\n", z.j);
                }
                else if (::_stricmp(key.c_str(), "K") == 0) {
                    if (!read_int_value(ifs, z.k, iline)) { error("failed reading K value"); return; }
                    info("    K = %d\n", z.k);
                }
                else if (::_stricmp(key.c_str(), "NODES") == 0) {
                    if (!read_int_value(ifs, z.nodes, iline)) { error("failed reading NODES value"); return; }
                    info("    NODES = %d\n", z.nodes);
                }
                else if (::_stricmp(key.c_str(), "ELEMENTS") == 0) {
                    if (!read_int_value(ifs, z.elements, iline)) { error("failed reading ELEMENTS value"); return; }
                    info("    ELEMENTS = %d\n", z.elements);
                }
                else if (::_stricmp(key.c_str(), "FACES") == 0) {
                    if (!read_int_value(ifs, z.faces, iline)) { error("failed reading FACES value"); return; }
                    info("    FACES = %d\n", z.faces);
                }
                else if (::_stricmp(key.c_str(), "TOTALNUMFACENODES") == 0) {
                    if (!read_int_value(ifs, z.num_face_nodes, iline)) { error("failed reading TotalNumFaceNodes value"); return; }
                    info("    TotalNumFaceNodes = %d\n", z.num_face_nodes);
                }
                else if (::_stricmp(key.c_str(), "TOTALNUMBOUNDARYCONNECTIONS") == 0) {
                    if (!read_int_value(ifs, z.num_bd_faces, iline)) { error("failed reading TotalNumBoundaryConnections value"); return; }
                    info("    TotalNumBoundaryConnections = %d\n", z.num_bd_faces);
                }
                // FACENEIGHBORMODE
                else if (::_stricmp(key.c_str(), "FACENEIGHBORMODE") == 0) {
                    // skip [...]
                    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '[');
                    if (ifs.peek() == '[')ifs.get();
                    ifs.ignore(std::numeric_limits<std::streamsize>::max(), ']');
                    if (ifs.peek() == ']')ifs.get();
                }
                // FACENEIGHBORCONNECTIONS
                else if (::_stricmp(key.c_str(), "FACENEIGHBORCONNECTIONS") == 0) {
                    int ival;
                    if (!read_int_value(ifs, ival, iline)) { error("failed reading FACENEIGHBORCONNECTIONS value"); return; }
                }
                // DT
                else if (::_stricmp(key.c_str(), "DT") == 0) {
                    // skip (xx,xxx,...)
                    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '(');
                    if (ifs.peek() == '(')ifs.get();
                    ifs.ignore(std::numeric_limits<std::streamsize>::max(), ')');
                    if (ifs.peek() == ')')ifs.get();
                }
                // DATAPACKING: BLOCK,POINT
                else if (::_stricmp(key.c_str(), "DATAPACKING") == 0) {
                    if (!read_string_value(ifs, value, iline)) { error("failed reading DATAPACKING value"); return; }
                    if (::_stricmp(value.c_str(), "BLOCK") == 0)z.is_block = true;
                    else if (::_stricmp(value.c_str(), "POINT") == 0)z.is_block = false;
                    else { error("invalid DATAPACKING: %s", value.c_str()); return; }
                }
                // VARLOCATION: ([]=xx,[]=xxx,...)
                else if (::_stricmp(key.c_str(), "VARLOCATION") == 0) {
                    if (!read_var_loc(ifs, z.var_loc, iline)) { error("failed reading VARLOCATION value"); return; }
                }
                // VARSHARELIST: ([]=xx,[]=xxx,...)
                else if (::_stricmp(key.c_str(), "VARSHARELIST") == 0) {
                    if (!read_var_share_list(ifs, z.var_share_list, zlast, iline)) { error("failed reading VARSHARELIST value"); return; }
                }
                // NV
                else if (::_stricmp(key.c_str(), "NV") == 0) {
                    int ival;
                    if (!read_int_value(ifs, ival, iline)) { error("failed reading NV value"); return; }
                }
                // CONNECTIVITYSHAREZONE
                else if (::_stricmp(key.c_str(), "CONNECTIVITYSHAREZONE") == 0) {
                    int ival;
                    if (!read_int_value(ifs, ival, iline)) { error("failed reading CONNECTIVITYSHAREZONE value"); return; }
                }
                // STRANDID
                else if (::_stricmp(key.c_str(), "STRANDID") == 0) {
                    int ival;
                    if (!read_int_value(ifs, ival, iline)) { error("failed reading STRANDID value"); return; }
                    info("    STRANDID = %d\n", ival);
                }
                // SOLUTIONTIME
                else if (::_stricmp(key.c_str(), "SOLUTIONTIME") == 0) {
                    if (!read_dbl_value(ifs, z.solution_time, iline)) { error("failed reading SOLUTIONTIME value"); return; }
                    info("    SolutionTime = %lf\n", z.solution_time);
                }
                // PARENTZONE
                else if (::_stricmp(key.c_str(), "PARENTZONE") == 0) {
                    int ival;
                    if (!read_int_value(ifs, ival, iline)) { error("failed reading PARENTZONE value"); return; }
                    info("    PARENTZONE = %d\n", ival);
                }
                // PASSIVEVARLIST: [...]
                else if (::_stricmp(key.c_str(), "PASSIVEVARLIST") == 0) {
                    // skip [...]
                    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '[');
                    if (ifs.peek() == '[')ifs.get();
                    ifs.ignore(std::numeric_limits<std::streamsize>::max(), ']');
                    if (ifs.peek() == ']')ifs.get();
                }
                // AUXDATA: NAME=STR
                else if (::_stricmp(key.c_str(), "AUXDATA") == 0) {
                    if (!read_key(ifs, key, iline)) { error("failed reading AUXDATA value"); return; }
                    if (!read_string_value(ifs, value, iline)) { error("failed reading AUXDATA value"); return; }
                    info("    AUXDATA %s = %s\n", key.c_str(), value.c_str());
                }
                else {
                    error("unknown key: %s\n", key.c_str()); return;
                }
            }

            // old format
            if (z.type == FEQUADRILATERAL) {
                if (z.i > 0)z.nodes    = z.i;
                if (z.j > 0)z.elements = z.j;
            }
            else if (z.type == ORDERED) {
                //if (z.i_ == 0)z.i_ = 1;
                //if (z.j_ == 0)z.j_ = 1;
                //if (z.k_ == 0)z.k_ = 1;
                int ni = std::max(z.i, 1);
                int nj = std::max(z.j, 1);
                int nk = std::max(z.k, 1);
                z.nodes    = ni * nj * nk;
                z.elements = std::max(ni - 1, 1) * std::max(nj - 1, 1) * std::max(nk - 1, 1);
            }

            // allocate
            z.vars.resize(variables_.size());
            for (size_t i = 0; i < variables_.size(); ++i) {
                if (z.var_loc[i] == NODAL)
                    z.vars[i].resize(z.nodes, 0);
                else
                    z.vars[i].resize(z.elements, 0);
            }

            // reading zone data

            skip_whitespace_and_comment(ifs, iline);

            info("    reading variables\n");

            // point
            if (!z.is_block) {
                for (int node = 0; node < z.nodes; ++node) {
                    for (size_t iv = 0; iv < variables_.size(); ++iv) {
                        ifs >> z.vars.at(iv).at(node);
                        if (!ifs.good())break;
                    }
                    if (!ifs.good())break;
                }
                if (!ifs.good()) { error("failed reading variable values"); return; }
            }
            // block
            else {
                for (auto& var : z.vars) {
                    skip_whitespace_and_comment(ifs, iline);
                    for (auto& x : var) {
                        ifs >> x;
                        if (!ifs.good())break;
                    }
                    if (!ifs.good())break;
                }
                if (!ifs.good()) { error("failed reading variable values"); return; }
            }

            // connection

            info("    reading connections\n");
            skip_whitespace_and_comment(ifs, iline);
            if      (z.type == FELINESEG) {
                // n11 n12
                // n21 n22
                // ...
                z.elem_nodes.reserve(z.elements, 2 * z.elements);
                for (int i = 0; i < z.elements; ++i) {
                    int nodes[2] = { 0 };
                    ifs >> nodes[0] >> nodes[1];
                    if (!ifs.good())break;

                    --nodes[0]; //? convert to zero-based index
                    --nodes[1];
                    z.elem_nodes.push_back(nodes);
                }
                if (!ifs.good()) { error("failed reading element-nodes data"); return; }
            }
            else if (z.type == FETRIANGLE) {
                // n11 n12 n13
                // n21 n22 n23
                // ...
                z.elem_nodes.reserve(z.elements, 3 * z.elements);
                for (int i = 0; i < z.elements; ++i) {
                    int nodes[3] = { 0 };
                    ifs >> nodes[0] >> nodes[1] >> nodes[2];
                    if (!ifs.good())break;

                    --nodes[0]; //? convert to zero-based index
                    --nodes[1];
                    --nodes[2];
                    z.elem_nodes.push_back(nodes);
                }
                if (!ifs.good()) { error("failed reading element-nodes data"); return; }
            }
            else if (z.type == FEQUADRILATERAL) {
                // n11 n12 n13 n14
                // n21 n22 n23 n24
                // ...
                z.elem_nodes.reserve(z.elements, 4 * z.elements);
                for (int i = 0; i < z.elements; ++i) {
                    int nodes[4] = { 0 };
                    ifs >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3];
                    if (!ifs.good())break;

                    --nodes[0]; //? convert to zero-based index
                    --nodes[1];
                    --nodes[2];
                    --nodes[3];
                    if (nodes[2] != nodes[3])
                        z.elem_nodes.push_back(nodes);
                    else
                        z.elem_nodes.push_back(3, nodes);
                }
                if (!ifs.good()) { error("failed reading element-nodes data"); return; }
            }
            else if (z.type == FETETRAHEDRON) {
                // n11 n12 n13 n14
                // n21 n22 n23 n24
                // ...
                z.elem_nodes.reserve(z.elements, 4 * z.elements);
                for (int i = 0; i < z.elements; ++i) {
                    int nodes[4] = { 0 };
                    ifs >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3];
                    if (!ifs.good())break;

                    --nodes[0]; //? convert to zero-based index
                    --nodes[1];
                    --nodes[2];
                    --nodes[3];
                    z.elem_nodes.push_back(nodes);
                }
                if (!ifs.good()) { error("failed reading element-nodes data"); return; }
            }
            else if (z.type == FEBRICK) {
                // n11 n12 n13 n14 n15 n16 n17 n18
                // n21 n22 n23 n24 n25 n26 n27 n28
                // ...
                z.elem_nodes.reserve(z.elements, 2 * z.elements);
                for (int i = 0; i < z.elements; ++i) {
                    int nodes[8] = { 0 };
                    ifs >> nodes[0] >> nodes[1] >> nodes[2] >> nodes[3] >> nodes[4] >> nodes[5] >> nodes[6] >> nodes[7];
                    if (!ifs.good())break;

                    --nodes[0]; //? convert to zero-based index
                    --nodes[1];
                    --nodes[2];
                    --nodes[3];
                    --nodes[4];
                    --nodes[5];
                    --nodes[6];
                    --nodes[7];
                    z.elem_nodes.push_back(nodes);
                }
                if (!ifs.good()) { error("failed reading element-nodes data"); return; }
            }
            else if (z.type == FEPOLYGON) {
                // face-nodes
                z.face_nodes.reserve(z.faces, 2 * z.faces);
                for (int i = 0; i < z.faces; ++i) {
                    int nodes[2] = { 0 };
                    ifs >> nodes[0] >> nodes[1];
                    if (!ifs.good())break;

                    --nodes[0]; //? convert to zero-based index
                    --nodes[1];
                    z.face_nodes.push_back(nodes);
                }
                if (!ifs.good()) { error("failed reading face-nodes data"); return; }

                // face-elems
                skip_whitespace_and_comment(ifs, iline);
                z.face_elems.resize(z.faces);
                for (auto& p : z.face_elems) {
                    ifs >> p.first >> p.second;
                    if (!ifs.good())break;
                    --p.first; --p.second;//? convert to zero-based index
                }
                if (!ifs.good()) { error("failed reading face-elements data"); return; }
            }
            else if (z.type == FEPOLYHEDRAL) {
                // node per face
                std::vector<int> nn2f(z.faces, 0);
                int len = 0;
                for (auto& n : nn2f) { ifs >> n; if (!ifs.good())break; len += n; }
                if (!ifs.good()) { error("failed reading face-nodes data"); return; }
                z.face_nodes.reserve(z.faces, len);

                skip_whitespace_and_comment(ifs, iline);

                // face-nodes
                for (int i = 0; i < z.faces; ++i) {
                    auto nodes = z.face_nodes.push_back(nn2f.at(i));
                    for (auto& n : nodes) {
                        ifs >> n;
                        --n;//? convert to zero-based index
                    }
                    if (!ifs.good())break;
                }
                if (!ifs.good()) { error("failed reading face-nodes data"); return; }

                // face-elems
                skip_whitespace_and_comment(ifs, iline);
                z.face_elems.resize(z.faces);
                for (auto& p : z.face_elems) {
                    ifs >> p.first >> p.second;
                    if (!ifs.good())break;
                    --p.first; --p.second;//? convert to zero-based index
                }
                if (!ifs.good()) { error("failed reading face-elements data"); return; }
            }

            info("  !!!OK!!!\n");

            //
            skip_whitespace_and_comment(ifs, iline);
            if (ifs.eof())break;

            if (!ifs.eof() && read_key(ifs, key, iline)) {
                if (::_stricmp(key.c_str(), "ZONE") != 0) {
                    error("unknown key: %s", key.c_str());
                }
            }

            ++zlast;
        } // next zone

        ifs.close();
        info("!!!OK!!!\n");
    }

}
