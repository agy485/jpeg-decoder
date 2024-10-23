// main.cpp
#include <stdio.h>
#include <jpeglib.h>
#include <setjmp.h>
#include <vector>
#include <map>
#include <algorithm>
#include <queue>
#include <climits>
#include <cmath>
#include <ostream>
#include <iostream>
#include <cstring>

struct my_error_mgr {
    struct jpeg_error_mgr pub;    // "public" fields
    jmp_buf setjmp_buffer;        // for return to caller
};

typedef struct my_error_mgr * my_error_ptr;

METHODDEF(void) my_error_exit(j_common_ptr cinfo)
{
    // cinfo->err 指向的是 jpeg_error_mgr
    my_error_ptr myerr = (my_error_ptr) cinfo->err;
    // 打印错误信息
    (*cinfo->err->output_message)(cinfo);
    // 跳转回 setjmp
    longjmp(myerr->setjmp_buffer, 1);
}

#define MAX_CLEN 16  // JPEG standard limit
#define MAX_SYMBOLS 257  // Maximum number of symbols

typedef struct {
    int symbol;
    long freq;
    int code_length;
    unsigned int code;
} Symbol;

typedef struct Node {
    long freq;
    int symbol; // Only valid for leaves
    struct Node *left;
    struct Node *right;
} Node;

// Comparator functions
int compare_symbols_by_freq(const void *a, const void *b);
int compare_symbols_by_freq_desc(const void *a, const void *b);
int compare_symbols_by_length_and_symbol(const void *a, const void *b);
int compare_nodes(const void *a, const void *b);

// Function prototypes
void assign_code_lengths(Node *node, int depth, Symbol symbols[], int num_symbols, int *max_code_length);
void free_huffman_tree(Node *node);

int find_common_prefix_length(const Symbol symbols[], int num_symbols) {
    if (num_symbols == 0) return 0;

    // 找到最大碼長
    int max_length = 17;
    for (int i = 0; i < num_symbols; i++) {
        if (symbols[i].code_length < max_length) {
            max_length = symbols[i].code_length;
        }
    }

    // 從最高位開始累加比較
    for (int bit_pos = 0; bit_pos < max_length; bit_pos++) {
        // 取第一個符號的當前位作為比較對象，並進行碼長檢查
        unsigned int prefix = 0;
        if (bit_pos < symbols[0].code_length) {
            prefix = (symbols[0].code >> (symbols[0].code_length - bit_pos - 1)) & 1;
        }

        bool mismatch = false;
        for (int i = 1; i < num_symbols; i++) {
            // 檢查當前符號碼長，避免比較超過符號有效位元長度的位
            unsigned int current_bit = 0;
            if (bit_pos < symbols[i].code_length) {
                current_bit = (symbols[i].code >> (symbols[i].code_length - bit_pos - 1)) & 1;
            }

            // 比較每個符號在當前位的值
            if (current_bit != prefix) {
                // 如果發現不一致，則返回當前的前綴長度
                printf("Mismatch found at bit position: %d\n", bit_pos);
                return bit_pos;
            }
        }

        // 在內層迴圈結束後，檢查所有符號移除前 bit_pos + 1 位後，剩餘的位元是否全為 1
        bool all_remaining_bits_are_one = false;
        for (int i = 0; i < num_symbols; i++) {
            int remaining_bits = symbols[i].code_length - (bit_pos + 1);
            if (remaining_bits <= 0) {
                // 沒有剩餘位元，繼續檢查下一個符號
                continue;
            }
            // 建立遮罩來取得剩餘位元
            unsigned int mask = (remaining_bits >= 32) ? ~0U : ((1U << remaining_bits) - 1);
            unsigned int remaining = symbols[i].code & mask;
            if (remaining == mask) {
                // 如果剩餘位元不是全 1，則不滿足條件
                all_remaining_bits_are_one = true;
                break;
            }
        }

        if (all_remaining_bits_are_one) {
            printf("All remaining bits are 1 after removing %d bits.\n", bit_pos + 1);
            return bit_pos;
        }
    }

    // 如果所有符號在所有位上都一致，則返回最大長度
    return max_length;
}


void remove_common_prefix(Symbol symbols[], int num_symbols, int common_prefix_length) {
    if (common_prefix_length == 0) return;

    // 將所有符號的編碼左移以去除共同前綴
    for (int i = 0; i < num_symbols; i++) {
        // symbols[i].code >>= common_prefix_length;
        symbols[i].code_length -= common_prefix_length;
    }
}

void generate_optimal_huffman_table(const long freq[257], JHUFF_TBL* htbl) {
    int num_symbols = 0;
    Symbol symbols[MAX_SYMBOLS];

    // Step 1: Collect symbols with non-zero frequency
    for (int i = 0; i < MAX_SYMBOLS; i++) {
        if (freq[i] > 0) {
            symbols[num_symbols].symbol = i;
            symbols[num_symbols].freq = freq[i];
            symbols[num_symbols].code_length = 0;
            symbols[num_symbols].code = 0;
            num_symbols++;
        }
    }

    if (num_symbols == 0) {
        // No symbols to encode
        memset(htbl->bits, 0, sizeof(htbl->bits));
        memset(htbl->huffval, 0, sizeof(htbl->huffval));
        return;
    }

    // Step 2: Sort symbols by increasing frequency
    qsort(symbols, num_symbols, sizeof(Symbol), compare_symbols_by_freq);

    // Step 3: Build the initial leaf nodes
    // Implement the standard Huffman algorithm

    // Create leaf nodes
    Node **nodes = (Node **)malloc(num_symbols * sizeof(Node *));
    for (int i = 0; i < num_symbols; i++) {
        nodes[i] = (Node *)malloc(sizeof(Node));
        nodes[i]->freq = symbols[i].freq;
        nodes[i]->symbol = symbols[i].symbol;
        nodes[i]->left = nodes[i]->right = NULL;
    }

    int num_nodes = num_symbols;

    // Build the Huffman tree
    while (num_nodes > 1) {
        // Sort nodes by frequency
        qsort(nodes, num_nodes, sizeof(Node *), compare_nodes);

        // Merge the two nodes with smallest frequencies
        Node *left = nodes[0];
        Node *right = nodes[1];

        Node *new_node = (Node *)malloc(sizeof(Node));
        new_node->freq = left->freq + right->freq;
        new_node->symbol = -1; // Not a leaf
        new_node->left = left;
        new_node->right = right;

        // Replace the two nodes with the new node
        nodes[0] = new_node;
        nodes[1] = nodes[num_nodes - 1];
        num_nodes--;
    }

    // Now assign code lengths
    int max_code_length = 0;
    assign_code_lengths(nodes[0], 0, symbols, num_symbols, &max_code_length);

    // Free the Huffman tree nodes
    free_huffman_tree(nodes[0]);
    free(nodes);

    // Adjust code lengths to fit within MAX_CLEN
    int bits[MAX_CLEN + 1] = {0}; // bits[1..16]

    for (int i = 0; i < num_symbols; i++) {
        int len = symbols[i].code_length;
        if (len > MAX_CLEN)
            len = MAX_CLEN;
        symbols[i].code_length = len;
        bits[len]++;
    }

    // Ensure that the total number of codes does not exceed the maximum allowed by the JPEG standard
    int total = 0;
    for (int i = MAX_CLEN; i > 0; i--) {
        total += bits[i];
        if (total > (1 << i)) {
            // Reduce bits[i] to ensure total codes <= 2^i
            int excess = total - (1 << i);
            bits[i] -= excess;
            bits[i - 1] += excess << 1; // Move excess codes to the next length
            total = 1 << i;
        }
    }

    // Re-assign code lengths based on the adjusted bits[]
    // First, sort symbols by frequency in descending order
    qsort(symbols, num_symbols, sizeof(Symbol), compare_symbols_by_freq_desc);

    // Now assign code lengths based on bits[]
    int idx = 0;
    for (int len = 1; len <= MAX_CLEN; len++) {
        int count = bits[len];
        for (int i = 0; i < count && idx < num_symbols; i++, idx++) {
            symbols[idx].code_length = len;
        }
    }

    // Sort symbols by code length and symbol value
    qsort(symbols, num_symbols, sizeof(Symbol), compare_symbols_by_length_and_symbol);

    // Generate the canonical Huffman codes
    unsigned int code = 0;
    int last_len = 0;
    for (int i = 0; i < num_symbols; i++) {
        int len = symbols[i].code_length;
        if (len > last_len) {
            code <<= (len - last_len);
            last_len = len;
        }
        symbols[i].code = code;
        code++;
    }

    int common_prefix_length = find_common_prefix_length(symbols, num_symbols);

    remove_common_prefix(symbols, num_symbols, common_prefix_length);

    // Build the Huffman table
    memset(htbl->bits, 0, sizeof(htbl->bits));
    memset(htbl->huffval, 0, sizeof(htbl->huffval));

    int hv_index = 0;
    for (int i = 0; i < num_symbols; i++) {
        int len = symbols[i].code_length;
        htbl->bits[len]++;
        htbl->huffval[hv_index++] = (unsigned char)symbols[i].symbol;
    }

    // Print out the code values for each symbol
    printf("Symbol\tFrequency\tCode Length\tCode\n");
    for (int i = 0; i < num_symbols; i++) {
        printf("%d\t%ld\t\t%d\t\t", symbols[i].symbol, symbols[i].freq, symbols[i].code_length);
        for (int j = symbols[i].code_length - 1; j >= 0; j--) {
            printf("%d", (symbols[i].code >> j) & 1);
        }
        printf("\n");
    }
}

// Comparator functions
int compare_symbols_by_freq(const void *a, const void *b) {
    const Symbol *sa = (const Symbol *)a;
    const Symbol *sb = (const Symbol *)b;
    if (sa->freq < sb->freq)
        return -1;
    else if (sa->freq > sb->freq)
        return 1;
    else
        return sa->symbol - sb->symbol;
}

int compare_symbols_by_freq_desc(const void *a, const void *b) {
    const Symbol *sa = (const Symbol *)a;
    const Symbol *sb = (const Symbol *)b;
    if (sa->freq > sb->freq)
        return -1;
    else if (sa->freq < sb->freq)
        return 1;
    else
        return sa->symbol - sb->symbol;
}

int compare_symbols_by_length_and_symbol(const void *a, const void *b) {
    const Symbol *sa = (const Symbol *)a;
    const Symbol *sb = (const Symbol *)b;
    if (sa->code_length != sb->code_length)
        return sa->code_length - sb->code_length;
    else
        return sa->symbol - sb->symbol;
}

int compare_nodes(const void *a, const void *b) {
    Node *na = *(Node **)a;
    Node *nb = *(Node **)b;
    if (na->freq < nb->freq)
        return -1;
    else if (na->freq > nb->freq)
        return 1;
    else
        return 0;
}

// Assign code lengths to the symbols
void assign_code_lengths(Node *node, int depth, Symbol symbols[], int num_symbols, int *max_code_length) {
    if (node->left == NULL && node->right == NULL) {
        // Leaf node
        for (int i = 0; i < num_symbols; i++) {
            if (symbols[i].symbol == node->symbol) {
                symbols[i].code_length = depth;
                if (depth > *max_code_length)
                    *max_code_length = depth;
                break;
            }
        }
    } else {
        if (node->left)
            assign_code_lengths(node->left, depth + 1, symbols, num_symbols, max_code_length);
        if (node->right)
            assign_code_lengths(node->right, depth + 1, symbols, num_symbols, max_code_length);
    }
}

// Free the Huffman tree nodes
void free_huffman_tree(Node *node) {
    if (node == NULL)
        return;
    free_huffman_tree(node->left);
    free_huffman_tree(node->right);
    free(node);
}

void print_huffman_table(const JHUFF_TBL* htbl) {
    std::cout << "Bits (碼長分布):" << std::endl;
    for (int i = 1; i <= 16; i++) {
        std::cout << "長度 " << i << ": " << (int)htbl->bits[i] << " 個符號" << std::endl;
    }

    std::cout << "Huffval (符號列表):" << std::endl;

    int total_symbols = 0;
    for (int i = 1; i <= 16; i++) {
        total_symbols += htbl->bits[i];
    }

    for (int i = 0; i < total_symbols; i++) {
        std::cout << "0x" << std::hex << (int)htbl->huffval[i] << " ";
    }
    std::cout << std::dec << std::endl;
}

void analyze_jpeg(const char* filename,
                  std::map<int, int>& dc_freq_Y,
                  std::map<std::pair<int, int>, int>& ac_freq_Y,
                  std::map<int, int>& dc_freq_CbCr,
                  std::map<std::pair<int, int>, int>& ac_freq_CbCr) {
    FILE* infile = fopen(filename, "rb");
    if (!infile) {
        fprintf(stderr, "无法打开文件 %s\n", filename);
        return;
    }

    struct jpeg_decompress_struct cinfo;
    struct my_error_mgr jerr;

    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;

    if (setjmp(jerr.setjmp_buffer)) {
        jpeg_destroy_decompress(&cinfo);
        fclose(infile);
        return;
    }

    // 初始化解压缩对象
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);

    // 读取 DCT 系数
    jvirt_barray_ptr* coef_arrays = jpeg_read_coefficients(&cinfo);

    // 初始化前一个块的 DC 系数
    std::vector<int> previous_dc(cinfo.num_components, 0);

    // 遍历每个组件（通常是 Y, Cb, Cr）
    for (int comp = 0; comp < cinfo.num_components; comp++) {
        jpeg_component_info* comp_info = cinfo.comp_info + comp;
        int width_in_blocks = comp_info->width_in_blocks;
        int height_in_blocks = comp_info->height_in_blocks;

        // 确定当前组件使用的频率表
        std::map<int, int>* dc_freq_ptr;
        std::map<std::pair<int, int>, int>* ac_freq_ptr;
        if (comp == 0) {
            // Y 组件
            dc_freq_ptr = &dc_freq_Y;
            ac_freq_ptr = &ac_freq_Y;
        } else {
            // Cb 和 Cr 组件
            dc_freq_ptr = &dc_freq_CbCr;
            ac_freq_ptr = &ac_freq_CbCr;
        }

        // 遍历每个块行
        for (int blk_y = 0; blk_y < height_in_blocks; blk_y++) {
            // 访问虚拟系数数组，获取当前行的系数
            JBLOCKARRAY coef_buffers = (cinfo.mem->access_virt_barray)
                ((j_common_ptr)&cinfo, coef_arrays[comp], blk_y, 1, FALSE);

            // 遍历当前行中的每个块
            for (int blk_x = 0; blk_x < width_in_blocks; blk_x++) {
                // 使用 JCOEF* 指针来引用 DCT 系数
                JCOEF* block = coef_buffers[0][blk_x];

                // 统计 DC 符号
                int dc_diff = block[0] - previous_dc[comp];
                previous_dc[comp] = block[0];
                int dc_category = (dc_diff == 0) ? 0 : (int)std::floor(std::log2(std::abs(dc_diff)) + 1);
                (*dc_freq_ptr)[dc_category]++;

                // 统计 AC 符号
                int run_length = 0;
                for (int i = 1; i < DCTSIZE2; i++) {
                    if (block[i] == 0) {
                        run_length++;
                        if (run_length == 16) {
                            // 插入 ZRL 符号
                            (*ac_freq_ptr)[{15, 0}]++;
                            run_length = 0;
                        }
                    } else {
                        int ac_size = (int)std::floor(std::log2(std::abs(block[i])) + 1);
                        // 插入 (Run, Size) 符号
                        (*ac_freq_ptr)[{run_length, ac_size}]++;
                        run_length = 0;
                    }
                }
                if (run_length > 0) {
                    // 插入 EOB 符号
                    (*ac_freq_ptr)[{0, 0}]++;
                }
            }
        }
    }

    // 输出 Y 组件的 DC 频率
    std::cout << "===== Y 组件 DC 符号频率 =====" << std::endl;
    for (const auto& pair : dc_freq_Y) {
        std::cout << "类别 " << pair.first << ": " << pair.second << " 次" << std::endl;
    }

    // 输出 Y 组件的 AC 频率
    std::cout << "===== Y 组件 AC 符号频率 =====" << std::endl;
    for (const auto& pair : ac_freq_Y) {
        int run = pair.first.first;
        int size = pair.first.second;
        int symbol;

        if (run == 0 && size == 0) {
            symbol = 0x00; // EOB
        } else if (run == 15 && size == 0) {
            symbol = 0xF0; // ZRL
        } else {
            symbol = (run << 4) | size;
        }

        std::cout << "符号 0x" << std::hex << symbol << std::dec
                  << " (Run=" << run << ", Size=" << size << "): "
                  << pair.second << " 次" << std::endl;
    }

    // 输出 CbCr 组件的 DC 频率
    std::cout << "===== CbCr 组件 DC 符号频率 =====" << std::endl;
    for (const auto& pair : dc_freq_CbCr) {
        std::cout << "类别 " << pair.first << ": " << pair.second << " 次" << std::endl;
    }

    // 输出 CbCr 组件的 AC 频率
    std::cout << "===== CbCr 组件 AC 符号频率 =====" << std::endl;
    for (const auto& pair : ac_freq_CbCr) {
        int run = pair.first.first;
        int size = pair.first.second;
        int symbol;

        if (run == 0 && size == 0) {
            symbol = 0x00; // EOB
        } else if (run == 15 && size == 0) {
            symbol = 0xF0; // ZRL
        } else {
            symbol = (run << 4) | size;
        }

        std::cout << "符号 0x" << std::hex << symbol << std::dec
                  << " (Run=" << run << ", Size=" << size << "): "
                  << pair.second << " 次" << std::endl;
    }

    // 清理
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
}

// 處理dc
JHUFF_TBL* create_huffman_table(j_compress_ptr cinfo, const std::map<int, int>& dc_freq_map) {
    // 分配新的 JHUFF_TBL 對象
    JHUFF_TBL* huff_tbl = jpeg_alloc_huff_table((j_common_ptr)cinfo);

    if (!huff_tbl) {
        fprintf(stderr, "無法分配霍夫曼表\n");
        return NULL;
    }

    // 初始化頻率數組
    long freq[257] = {0}; // 索引 0 - 256

    // 將符號頻率填入 freq 數組
    for (const auto& pair : dc_freq_map) {
        int symbol = pair.first;
        int frequency = pair.second;

        // 確保符號值在有效範圍内（DC 分類通常為 0 - 11）
        if (symbol < 0 || symbol > 11) {
            fprintf(stderr, "DC 符號超出範圍: %d\n", symbol);
            continue;
        }

        freq[symbol] = frequency;
    }

    generate_optimal_huffman_table(freq, huff_tbl);

    return huff_tbl;
}

// 處理ac
JHUFF_TBL* create_huffman_table_ac(j_compress_ptr cinfo, const std::map<std::pair<int, int>, int>& ac_freq_map) {
    JHUFF_TBL* huff_tbl = jpeg_alloc_huff_table((j_common_ptr)cinfo);

    if (!huff_tbl) {
        fprintf(stderr, "无法分配霍夫曼表\n");
        return NULL;
    }

    // 初始化频率數組
    long freq[257] = {0}; // 索引 0 - 256

    // 將 (Run, Size) 符號映射到符號值，並填入 freq 數組
    for (const auto& pair : ac_freq_map) {
        int run = pair.first.first;
        int size = pair.first.second;
        int symbol;

        if (run == 0 && size == 0) {
            symbol = 0x00; // EOB
        } else if (run == 15 && size == 0) {
            symbol = 0xF0; // ZRL
        } else {
            symbol = (run << 4) | size;
        }

        // 確保符號值在 0 - 255 之間
        if (symbol < 0 || symbol > 255) {
            fprintf(stderr, "AC 符號超出範圍: run=%d, size=%d, symbol=%d\n", run, size, symbol);
            continue;
        }

        freq[symbol] = pair.second;
    }

    // 确保所有可能的 AC 符号频率至少为 1
    for (int run = 0; run <= 15; run++) {
        for (int size = 1; size <= 10; size++) {
            int symbol = (run << 4) | size;
            if (freq[symbol] == 0) {
                freq[symbol] = 1;
            }
        }
    }
    // 对特殊符号 EOB 和 ZRL 也进行处理
    if (freq[0x00] == 0) {
        freq[0x00] = 1; // EOB
    }
    if (freq[0xF0] == 0) {
        freq[0xF0] = 1; // ZRL
    }

    generate_optimal_huffman_table(freq, huff_tbl);

    return huff_tbl;
}

void write_jpeg_with_custom_huffman(const char* input_filename, const char* output_filename,
                                    const std::map<int, int>& dc_freq_Y, const std::map<std::pair<int, int>, int>& ac_freq_Y,
                                    const std::map<int, int>& dc_freq_CbCr, const std::map<std::pair<int, int>, int>& ac_freq_CbCr) {
    FILE* infile = fopen(input_filename, "rb");
    if (!infile) {
        fprintf(stderr, "無法打開文件 %s\n", input_filename);
        return;
    }

    FILE* outfile = fopen(output_filename, "wb");
    if (!outfile) {
        fprintf(stderr, "無法創建文件 %s\n", output_filename);
        fclose(infile);
        return;
    }

    struct jpeg_decompress_struct cinfo;
    struct jpeg_compress_struct out_cinfo;
    struct my_error_mgr jerr;

    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;

    out_cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;

    if (setjmp(jerr.setjmp_buffer)) {
        jpeg_destroy_decompress(&cinfo);
        jpeg_destroy_compress(&out_cinfo);
        fclose(infile);
        fclose(outfile);
        return;
    }

    // 初始化解壓縮對象
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);

    // 開始解壓縮
    jpeg_start_decompress(&cinfo);

    // 初始化壓縮對象
    jpeg_create_compress(&out_cinfo);
    jpeg_stdio_dest(&out_cinfo, outfile);

    // 設置壓縮參數
    out_cinfo.image_width = cinfo.image_width;
    out_cinfo.image_height = cinfo.image_height;
    out_cinfo.input_components = 3;  // JPEG 預設 YCbCr 色彩空間有三個成分（Y, Cb, Cr）
    // out_cinfo.in_color_space = JCS_YCbCr;
    out_cinfo.in_color_space = JCS_RGB;
    // cinfo.out_color_space = out_cinfo.in_color_space;

    jpeg_set_defaults(&out_cinfo);

    // 禁用優化編碼，以使用自定義霍夫曼表
    out_cinfo.optimize_coding = FALSE;

 // 为 Y 组件生成自定义霍夫曼表
    JHUFF_TBL* custom_dc_tbl_Y = create_huffman_table(&out_cinfo, dc_freq_Y);
    JHUFF_TBL* custom_ac_tbl_Y = create_huffman_table_ac(&out_cinfo, ac_freq_Y);

    // 为 CbCr 组件生成自定义霍夫曼表
    JHUFF_TBL* custom_dc_tbl_C = create_huffman_table(&out_cinfo, dc_freq_CbCr);
    JHUFF_TBL* custom_ac_tbl_C = create_huffman_table_ac(&out_cinfo, ac_freq_CbCr);

    // 检查是否成功
    if (!custom_dc_tbl_Y || !custom_ac_tbl_Y || !custom_dc_tbl_C || !custom_ac_tbl_C) {
        fprintf(stderr, "自定义霍夫曼表生成失败\n");
        jpeg_destroy_decompress(&cinfo);
        jpeg_destroy_compress(&out_cinfo);
        fclose(infile);
        fclose(outfile);
        return;
    }

    // 分配并设置自定义 DC 霍夫曼表
    out_cinfo.dc_huff_tbl_ptrs[0] = jpeg_alloc_huff_table((j_common_ptr)&out_cinfo);
    memcpy(out_cinfo.dc_huff_tbl_ptrs[0]->bits, custom_dc_tbl_Y->bits, sizeof(custom_dc_tbl_Y->bits));
    memcpy(out_cinfo.dc_huff_tbl_ptrs[0]->huffval, custom_dc_tbl_Y->huffval, sizeof(custom_dc_tbl_Y->huffval));

    out_cinfo.dc_huff_tbl_ptrs[1] = jpeg_alloc_huff_table((j_common_ptr)&out_cinfo);
    memcpy(out_cinfo.dc_huff_tbl_ptrs[1]->bits, custom_dc_tbl_C->bits, sizeof(custom_dc_tbl_C->bits));
    memcpy(out_cinfo.dc_huff_tbl_ptrs[1]->huffval, custom_dc_tbl_C->huffval, sizeof(custom_dc_tbl_C->huffval));

    // 分配并设置自定义 AC 霍夫曼表
    out_cinfo.ac_huff_tbl_ptrs[0] = jpeg_alloc_huff_table((j_common_ptr)&out_cinfo);
    memcpy(out_cinfo.ac_huff_tbl_ptrs[0]->bits, custom_ac_tbl_Y->bits, sizeof(custom_ac_tbl_Y->bits));
    memcpy(out_cinfo.ac_huff_tbl_ptrs[0]->huffval, custom_ac_tbl_Y->huffval, sizeof(custom_ac_tbl_Y->huffval));

    out_cinfo.ac_huff_tbl_ptrs[1] = jpeg_alloc_huff_table((j_common_ptr)&out_cinfo);
    memcpy(out_cinfo.ac_huff_tbl_ptrs[1]->bits, custom_ac_tbl_C->bits, sizeof(custom_ac_tbl_C->bits));
    memcpy(out_cinfo.ac_huff_tbl_ptrs[1]->huffval, custom_ac_tbl_C->huffval, sizeof(custom_ac_tbl_C->huffval));

    // 对组件分别设置霍夫曼表编号
    out_cinfo.comp_info[0].dc_tbl_no = 0; // Y 组件使用表 0
    out_cinfo.comp_info[0].ac_tbl_no = 0;

    out_cinfo.comp_info[1].dc_tbl_no = 1; // Cb 组件使用表 1
    out_cinfo.comp_info[1].ac_tbl_no = 1;

    out_cinfo.comp_info[2].dc_tbl_no = 1; // Cr 组件使用表 1
    out_cinfo.comp_info[2].ac_tbl_no = 1;

    jpeg_start_compress(&out_cinfo, TRUE);

    JSAMPARRAY buffer = (*cinfo.mem->alloc_sarray)
        ((j_common_ptr)&cinfo, JPOOL_IMAGE, cinfo.output_width * cinfo.output_components, 1);
    while (cinfo.output_scanline < cinfo.output_height) {
        jpeg_read_scanlines(&cinfo, buffer, 1);
        jpeg_write_scanlines(&out_cinfo, buffer, 1);
    }

    // 完成壓縮和解壓縮
    jpeg_finish_compress(&out_cinfo);
    jpeg_destroy_compress(&out_cinfo);

    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);

    // 获取 infile 和 outfile 的文件大小
    fseek(infile, 0, SEEK_END);
    long infile_size = ftell(infile);
    fseek(outfile, 0, SEEK_END);
    long outfile_size = ftell(outfile);

    printf("输入文件大小: %ld 字节\n", infile_size);
    printf("输出文件大小: %ld 字节\n", outfile_size);

    fclose(infile);
    fclose(outfile);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "使用方法: %s <輸入JPEG> <輸出JPEG>\n", argv[0]);
        return 1;
    }

    const char* input_filename = argv[1];
    const char* output_filename = argv[2];

    // 步驟1：分析JPEG並統計符號頻率
    std::map<int, int> dc_freq_Y;
    std::map<std::pair<int, int>, int> ac_freq_Y;
    std::map<int, int> dc_freq_CbCr;
    std::map<std::pair<int, int>, int> ac_freq_CbCr;

    // 分析 JPEG 图像
    analyze_jpeg(input_filename, dc_freq_Y, ac_freq_Y, dc_freq_CbCr, ac_freq_CbCr);

    // 步驟2和步驟3：生成自定義霍夫曼表並重新編碼保存JPEG
    write_jpeg_with_custom_huffman(input_filename, output_filename,
                                   dc_freq_Y, ac_freq_Y, dc_freq_CbCr, ac_freq_CbCr);

    return 0;
}
