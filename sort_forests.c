#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/sendfile.h>

#define _FILE_OFFSET_BITS 64

#include "utils.h"
#include "sglib.h"
#include "progressbar.h"

#ifndef ADD_DIFF_TIME
#define ADD_DIFF_TIME(tbegin, tend)  ((tend.tv_sec - tstart.tv_sec) + 1e-6*(tend.tv_usec - tstart.tv_usec))
#endif

/* Make sure LOCATIONS_FILENAME is a multiple of 8 */
#define LOCATIONS_FILENAME_SIZE  (32)
struct locations{
    int64_t forestid; /* forest id for the tree */
    int64_t tree_root; /* tree root id (id for z=0 halo) */
    int64_t fileid;/* fileid, 1d index for the 3-d decomposition by consistent-trees */
    int64_t offset;/* the offset in tree_* file. */
    int64_t bytes;/* number of bytes in the file. MUST BE SIGNED type */
    char filename[LOCATIONS_FILENAME_SIZE];/* filename where the tree is written */
};

void usage(int argc, char **argv)
{
    (void) argc;
    fprintf(stderr,"USAGE: %s <input consistent-trees directory> <output_consistent_trees_directory> \n",
            argv[0]);
}    

int64_t read_forests(const char *filename, int64_t **f, int64_t **t)
{
    char buffer[MAXBUFSIZE];
    const char comment = '#';
    /* By passing the comment character, getnumlines
       will return the actual number of lines, ignoring
       the first header line. 
     */

    const int64_t ntrees = getnumlines(filename, comment);
    *f = my_malloc(sizeof(int64_t), ntrees);
    *t = my_malloc(sizeof(int64_t), ntrees);

    int64_t *forests    = *f;
    int64_t *tree_roots = *t;
    
    int64_t ntrees_found = 0;
    FILE *fp = my_fopen(filename, "r");
    while(fgets(buffer, MAXBUFSIZE, fp) != NULL) {
        if(buffer[0] == comment) {
            continue;
        } else {
            const int nitems_expected = 2;
            XASSERT(ntrees_found < ntrees,
                    "ntrees=%"PRId64" should be less than ntrees_found=%"PRId64"\n", ntrees, ntrees_found);            
            int nitems = sscanf(buffer, "%"SCNd64" %"SCNd64, &(tree_roots[ntrees_found]), &(forests[ntrees_found]));
            XASSERT(nitems == nitems_expected,
                    "Expected to parse %d long integers but found `%s' in the buffer. nitems = %d \n",
                    nitems_expected, buffer, nitems);
            ntrees_found++;
        }
    }
    XASSERT(ntrees == ntrees_found,
            "ntrees=%"PRId64" should be equal to ntrees_found=%"PRId64"\n", ntrees, ntrees_found);
    fclose(fp);

    return ntrees_found;
}    


int64_t read_locations(const char *filename, const int64_t ntrees, struct locations *l, int *num_files, int *box_div)
{
    char buffer[MAXBUFSIZE];
    int max_fileid = 0;
    const char comment = '#';
    /* By passing the comment character, getnumlines
       will return the actual number of lines, ignoring
       the first header line. 
     */

    struct locations *locations = l;

    int64_t ntrees_found = 0;
    FILE *fp = my_fopen(filename, "r");
    while(fgets(buffer, MAXBUFSIZE, fp) != NULL) {
        if(buffer[0] == comment) {
            continue;
        } else {
            const int nitems_expected = 4;
            char linebuf[MAXLEN];
            XASSERT(ntrees_found < ntrees,
                    "ntrees=%"PRId64" should be less than ntrees_found=%"PRId64"\n",
                    ntrees, ntrees_found);            
            int nitems = sscanf(buffer, "%"SCNd64" %"SCNd64 " %"SCNd64 "%s", &(locations[ntrees_found].tree_root),
                                &(locations[ntrees_found].fileid), &(locations[ntrees_found].offset), linebuf);

            /* The filename is separated out to save memory but I want to ensure that the actual filename does
               not get truncated. The filename field might actually be removed later. */
            my_snprintf(locations[ntrees_found].filename, LOCATIONS_FILENAME_SIZE, "%s", linebuf);
            XASSERT(nitems == nitems_expected,
                    "Expected to parse two long integers but found `%s' in the buffer\n",
                    buffer);
            ntrees_found++;
        }
    }
    XASSERT(ntrees == ntrees_found, "ntrees=%"PRId64" should be equal to ntrees_found=%"PRId64"\n", ntrees, ntrees_found);
    fclose(fp);

    for(int i=0;i<ntrees_found;i++){
        if (locations[i].fileid > max_fileid) {
            max_fileid = locations[i].fileid;
        }
    }

    /* number of files is one greater from 0 based indexing of C files */
    *num_files = max_fileid + 1;
    const int box_divisions = (int) round(cbrt(*num_files));
    const int box_cube = box_divisions * box_divisions * box_divisions;
    XASSERT( (box_cube) == (*num_files),
             "box_divisions^3=%d should be equal to nfiles=%d\n",
             box_cube, *num_files);
    *box_div = box_divisions;
    return ntrees_found;
}    


void sort_forests(const int64_t ntrees, int64_t *forests, int64_t *tree_roots)
{

#define MULTIPLE_ARRAY_EXCHANGER(type,a,i,j) { SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64_t, tree_roots,i,j); \
        SGLIB_ARRAY_ELEMENTS_EXCHANGER(int64_t, forests, i, j) }
    
    SGLIB_ARRAY_QUICK_SORT(int64_t, tree_roots, ntrees, SGLIB_NUMERIC_COMPARATOR , MULTIPLE_ARRAY_EXCHANGER);
#undef MULTIPLE_ARRAY_EXCHANGER

}    

int compare_locations_tree_roots(const void *l1, const void *l2)
{
    const struct locations *aa = (const struct locations *) l1;
    const struct locations *bb = (const struct locations *) l2;
    return (aa->tree_root < bb->tree_root) ? -1:1;
}

int compare_locations_fid(const void *l1, const void *l2)
{
    const struct locations *aa = (const struct locations *) l1;
    const struct locations *bb = (const struct locations *) l2;
    return (aa->forestid < bb->forestid) ? -1:1;
}

int compare_locations_file_offset(const void *l1, const void *l2)
{
    const struct locations *aa = (const struct locations *) l1;
    const struct locations *bb = (const struct locations *) l2;

    const int file_id_cmp = (aa->fileid == bb->fileid) ? 0:((aa->fileid < bb->fileid) ? -1:1);
    if(file_id_cmp == 0) {
        /* trees are in same file -> sort by offset */
        return (aa->offset < bb->offset) ? -1:1;
    } else {
        return file_id_cmp;
    }
        
    return 0;
}

int compare_locations_fid_file_offset(const void *l1, const void *l2)
{
    const struct locations *aa = (const struct locations *) l1;
    const struct locations *bb = (const struct locations *) l2;

    if(aa->forestid != bb->forestid) {
        return (aa->forestid < bb->forestid) ? -1:1;
    } else {
        /* The trees are in the same forest. Check filename */
        const int file_id_cmp = (aa->fileid == bb->fileid) ? 0:((aa->fileid < bb->fileid) ? -1:1);
        if(file_id_cmp == 0) {
            /* trees are in same file -> sort by offset */
            return (aa->offset < bb->offset) ? -1:1;
        } else {
            return file_id_cmp;
        }
    }
        
    return 0;
}

void sort_locations_on_treeroot(const int64_t ntrees, struct locations *locations)
{
    qsort(locations, ntrees, sizeof(*locations), compare_locations_tree_roots);
}    

void sort_locations_file_offset(const int64_t ntrees, struct locations *locations)
{
    qsort(locations, ntrees, sizeof(*locations), compare_locations_file_offset);
}

void sort_locations_on_fid(const int64_t ntrees, struct locations *locations)
{
    qsort(locations, ntrees, sizeof(*locations), compare_locations_fid);
}    

void sort_locations_on_fid_file_offset(const int64_t ntrees, struct locations *locations)
{
    qsort(locations, ntrees, sizeof(*locations), compare_locations_fid_file_offset);
}    


void assign_forest_ids(const int64_t ntrees, struct locations *locations, int64_t *forests, int64_t *tree_roots)
{
    /* Sort forests by tree roots -> necessary for assigning forest ids */
    sort_forests(ntrees, forests, tree_roots);
    sort_locations_on_treeroot(ntrees, locations);    
    
    /* forests and tree_roots are sorted together, on tree_roots */
    /* locations is sorted on tree roots */
    for(int64_t i=0;i<ntrees;i++) {
        XASSERT(tree_roots[i] == locations[i].tree_root,
                "tree roots[%"PRId64"] = %"PRId64" does not equal tree roots in locations = %"PRId64"\n",
                i, tree_roots[i], locations[i].tree_root);
        locations[i].forestid = forests[i];
    }
}    

void assign_trees_in_forest_to_same_file(const int64_t ntrees, struct locations *locations, struct locations *new_locations, const int nfiles, const int BOX_DIVISIONS)
{
    sort_locations_on_fid_file_offset(ntrees, locations);
    sort_locations_on_fid_file_offset(ntrees, new_locations);
    
    int64_t *histogram_fileids = my_calloc(sizeof(*histogram_fileids), nfiles);

    /* the fun begins here -> in case, assign all trees from a forest into the same file */
    int64_t start_forestid = locations[0].forestid;
    int64_t min_fileid = locations[0].fileid;
    int64_t max_fileid = locations[0].fileid;
    int64_t start_index_forest = 0;
    int64_t end_index_forest = 1;
    int64_t num_trees_moved = 0;

    fprintf(stderr, ANSI_COLOR_MAGENTA"Assigning all trees in a forest into the same file...."ANSI_COLOR_RESET"\n");
    /* setup the progressbar */
    int interrupted=0;
    init_my_progressbar(ntrees, &interrupted);
   
    for(int64_t i=1;i<ntrees;i++) {
        my_progressbar(i, &interrupted);
        if(locations[i].forestid == start_forestid) {
            if(locations[i].fileid < min_fileid) {
                min_fileid = locations[i].fileid;
            }
            if(locations[i].fileid > min_fileid) {
                max_fileid = locations[i].fileid;
            }
            end_index_forest++;
            continue;
        } else {
            if(min_fileid == max_fileid) {
                for(int64_t j=start_index_forest;j<end_index_forest;j++) {
                    new_locations[j].fileid = min_fileid;
                }
            } else {
                /* fprintf(stderr,"For forest id = %"PRId64" trees are stored in separate files (min, max) = (%"PRId64", %"PRId64")\n", */
                /*         start_forestid, min_fileid, max_fileid); */
                /* interrupted=1; */

                /* create a histogram of the fileids */
                memset(histogram_fileids, 0, sizeof(*histogram_fileids) * nfiles);
                for(int64_t j=start_index_forest;j<end_index_forest;j++) {
                    histogram_fileids[locations[j].fileid]++;
                }

                int64_t max_common_value = 0;
                int max_common_fileid = 0;
                for(int j=0;j<nfiles;j++) {
                    if(histogram_fileids[j] > max_common_value) {
                        max_common_value = histogram_fileids[j];
                        max_common_fileid = j;
                    }
                }

                const int ii = max_common_fileid/(BOX_DIVISIONS*BOX_DIVISIONS);
                const int jj = (max_common_fileid%((int64_t)(BOX_DIVISIONS*BOX_DIVISIONS)))/BOX_DIVISIONS;
                const int kk = max_common_fileid%((int64_t)BOX_DIVISIONS);
                for(int64_t j=start_index_forest;j<end_index_forest;j++) {
                    new_locations[j].fileid = max_common_fileid;
                    my_snprintf(new_locations[j].filename, LOCATIONS_FILENAME_SIZE, "tree_%d_%d_%d.dat", ii, jj, kk);
                    if(new_locations[j].fileid != locations[j].fileid) {
                        num_trees_moved++;
                        /* fprintf(stderr,"Moved tree = %10"PRId64" from fileid=%3"PRId64" to fileid=%3"PRId64"\n",locations[j].tree_root, locations[j].fileid, new_locations[j].fileid); */
                        /* interrupted=1; */
                    }
                }
            }

            start_forestid = locations[i].forestid;
            start_index_forest = i;
            end_index_forest   = i+1;
            min_fileid = locations[i].fileid;
            max_fileid = locations[i].fileid;
        }
    }
    finish_myprogressbar(&interrupted);
    free(histogram_fileids);
    if(num_trees_moved > 0) {
        fprintf(stderr,"Number of trees moved into different files = %"PRId64"\n",num_trees_moved);
    }
    fprintf(stderr, ANSI_COLOR_GREEN"Assigning all trees in a forest into the same file.......done"ANSI_COLOR_RESET"\n\n");
}


int64_t compute_numbytes_with_off(const int64_t off, const int64_t start)
{
    return off - start; /* Or should there be a -1?*/
}


int64_t compute_numbytes(FILE *fp, const int64_t start)
{
    const int64_t off = ftello(fp);
    return compute_numbytes_with_off(off, start);
}

int64_t write_forests_and_locations(const char *filename, const int64_t ntrees, const struct locations *locations)
{
    int64_t bytes = 0;
    FILE *fp = my_fopen(filename,"w");
    bytes += fprintf(fp,"############## Value-Added Consistent-Trees locations.dat file. #####################\n");
    bytes += fprintf(fp,"##    ForestID     TreeRootID     FileID    Offset        Filename          Bytes    \n");
    bytes += fprintf(fp,"#####################################################################################\n");
    
    for(int64_t i=0;i<ntrees;i++) {
        bytes += fprintf(fp," %13"PRId64" %13"PRId64" %6"PRId64" %13"PRId64"  %s %16"PRId64"\n",
                         locations[i].forestid, locations[i].tree_root, locations[i].fileid,
                         locations[i].offset, locations[i].filename, locations[i].bytes);
    }
    fclose(fp);
    return bytes;
}


int main(int argc, char **argv)
{
    char *input_dir, *output_dir;
    if(argc != 3) {
        usage(argc, argv);
        return EXIT_FAILURE;
    } else {
        input_dir  = argv[1];
        output_dir = argv[2];
    }
    
    if(strcmp(input_dir, output_dir) == 0) {
        fprintf(stderr,"ERROR: Input and output directories are the same..exiting\n");
        return EXIT_FAILURE;
    }

    struct timeval tstart, tend;
    gettimeofday(&tstart, NULL);
    char locations_filename[MAXLEN], forests_filename[MAXLEN];
    int64_t *forests=NULL, *tree_roots=NULL;
    my_snprintf(locations_filename, MAXLEN, "%s/locations.dat", input_dir);
    my_snprintf(forests_filename, MAXLEN, "%s/forests.list", input_dir);
    fprintf(stderr, ANSI_COLOR_MAGENTA"Reading forests...."ANSI_COLOR_RESET"\n");
    const int64_t ntrees = read_forests(forests_filename, &forests, &tree_roots);
    fprintf(stderr, ANSI_COLOR_GREEN"Reading forests......done"ANSI_COLOR_RESET"\n\n");
    /* fprintf(stderr, "Number of trees = %"PRId64"\n\n",ntrees); */

    struct locations *locations = my_malloc(sizeof(*locations), ntrees);
    int nfiles = 0, BOX_DIVISIONS=0;
    fprintf(stderr, ANSI_COLOR_MAGENTA"Reading locations...."ANSI_COLOR_RESET"\n");
    const int64_t ntrees_loc = read_locations(locations_filename, ntrees, locations, &nfiles, &BOX_DIVISIONS);
    fprintf(stderr, ANSI_COLOR_GREEN"Reading locations......done"ANSI_COLOR_RESET"\n\n");
    XASSERT(ntrees == ntrees_loc,
            "ntrees=%"PRId64" should be equal to ntrees_loc=%"PRId64"\n",
            ntrees, ntrees_loc);    

    /* the following function will sort locations and forests based on tree root id*/
    assign_forest_ids(ntrees, locations, forests, tree_roots);

    /* Forests are now contained inside locations -> free the pointers */
    free(forests);free(tree_roots);


    FILE **tree_outputs = my_malloc(sizeof(FILE *), nfiles);
    FILE **tree_inputs  = my_malloc(sizeof(FILE *), nfiles);

    int *tree_inputs_fd = my_malloc(sizeof(*tree_inputs_fd), nfiles);
    int *tree_outputs_fd = my_malloc(sizeof(*tree_outputs_fd), nfiles);    

    XASSERT(sizeof(off_t) == 8,
            "File offset bits must be 64\n"
            "Please ensure "ANSI_COLOR_RED"#define _FILE_OFFSET_BITS 64"ANSI_COLOR_RESET" is present\n");

    off_t *tree_outputs_fd_offset = my_malloc(sizeof(*tree_outputs_fd_offset), nfiles);


    int64_t *tree_counts = my_calloc(sizeof(*tree_counts), nfiles);
    int64_t *inp_file_sizes = my_calloc(sizeof(*inp_file_sizes), nfiles);
    char buffer[MAXLEN];
    for (int i=0; i<BOX_DIVISIONS; i++) {
        for (int j=0; j<BOX_DIVISIONS; j++) {
            for(int k=0; k<BOX_DIVISIONS; k++) {
                my_snprintf(buffer,MAXLEN,"%s/tree_%d_%d_%d.dat", input_dir, i, j, k);
                int id = id = i*BOX_DIVISIONS*BOX_DIVISIONS + j*BOX_DIVISIONS + k;
                tree_inputs[id]  = my_fopen(buffer, "r");
                
                XASSERT(setvbuf(tree_inputs[id], NULL, _IONBF, 0) == 0,
                        "Could not set unbuffered fgets");
                my_fseek(tree_inputs[id],0L, SEEK_END);
                inp_file_sizes[id] = ftello(tree_inputs[id]);
                rewind(tree_inputs[id]);

                tree_inputs_fd[id]  = fileno(tree_inputs[id]);

                my_snprintf(buffer,MAXLEN,"%s/tree_%d_%d_%d.dat", output_dir, i, j, k);
                unlink(buffer);
                tree_outputs[id] = my_fopen(buffer, "w");
                /* setbuf(tree_outputs[id], _IOFBF); */
                tree_outputs_fd[id] = fileno(tree_outputs[id]);
            }
        }
    }


    /* the following function will sort locations based on 1) filename 2) offsets */
    sort_locations_file_offset(ntrees, locations);

    /* holder to check later that bytes have been assigned */
    for(int64_t i=0;i<ntrees;i++) {
        locations[i].bytes = -1;/* Make sure bytes is a signed type! */
    }

    /* Create a copy of current locations */    
    struct locations *new_locations = my_malloc(sizeof(*new_locations), ntrees);
    assert(sizeof(*new_locations) == sizeof(*locations) && "locations struct is varying in size! The sky is falling!!");
    memcpy(new_locations, locations, sizeof(*locations) * ntrees);

    /* figure out the byte size for each tree */
    int64_t start = locations[0].offset;
    int64_t start_fileid = locations[0].fileid;

    /* tree_roots are 64 bit integers -> max digits in decimal = log10(2^64) < 20.
       Add 1 char for +-, in case consistent tree changes. and then strlen('#tree ')
       and the previous \n. I need to read up to previous newline.
    */
    const int64_t guess_max_linesize = 20 + 1 + 6 + 1;
    fprintf(stderr, ANSI_COLOR_MAGENTA"Calculating the number of bytes for each tree...."ANSI_COLOR_RESET"\n");
    /* setup the progressbar */
    int interrupted=0;
    init_my_progressbar(ntrees, &interrupted);

    for(int64_t i=1;i<=ntrees-1;i++) {
        my_progressbar(i, &interrupted);
        const int64_t fileid = locations[i].fileid;
        
        /* Are we starting on a new file ?*/
        if(start_fileid != fileid) {
            /* fill out the bytes for the last tree in the previous file */
            const int64_t num_bytes = compute_numbytes_with_off(inp_file_sizes[start_fileid], start);
            locations[i-1].bytes = num_bytes;
            new_locations[i-1].bytes = num_bytes;

            /* now we reset the start fields */
            start = locations[i].offset;
            start_fileid = locations[i].fileid;
            continue;
        }
        const int64_t current_offset_guess = locations[i].offset - guess_max_linesize;
        my_fseek(tree_inputs[fileid], current_offset_guess, SEEK_SET);
        while(1) {
            const int a = fgetc(tree_inputs[fileid]);
            if(a == EOF) {
                fprintf(stderr,"Encountered EOF while looking for end of current tree\n");
                exit(EXIT_FAILURE);
            }
            const unsigned char c = (unsigned char) a;
            if(c == '\n') {
                const int64_t num_bytes = compute_numbytes(tree_inputs[start_fileid], start);
                locations[i-1].bytes = num_bytes;
                new_locations[i-1].bytes = num_bytes;
                /* fprintf(stderr,"%"PRId64"\n",num_bytes); */
                start = locations[i].offset;
                break;
            }
        }
    }

    /* fill out the bytes for the last tree */
    {
        start = locations[ntrees-1].offset;
        const int64_t fileid = locations[ntrees-1].fileid;
        my_fseek(tree_inputs[fileid], 0L, SEEK_END);
        const int64_t num_bytes = compute_numbytes(tree_inputs[fileid], start);
        locations[ntrees-1].bytes = num_bytes;
        new_locations[ntrees-1].bytes = num_bytes;
    }
    finish_myprogressbar(&interrupted);        
    fprintf(stderr, ANSI_COLOR_GREEN"Calculating the number of bytes for each tree.....done"ANSI_COLOR_RESET"\n\n");

    for(int64_t i=ntrees-1;i>=0;i--) {
        XASSERT(locations[i].bytes > 0,
                "locations[%"PRId64"].bytes = %"PRId64" should be positive\n",
                i,locations[i].bytes);

        XASSERT(new_locations[i].bytes == locations[i].bytes,
                "locations[%"PRId64"].bytes = %"PRId64" should be equal new_locations->bytes = %"PRId64"\n",
                i,locations[i].bytes,new_locations[i].bytes);
        XASSERT(strncmp(new_locations[i].filename, locations[i].filename, LOCATIONS_FILENAME_SIZE) == 0,
                "new_locations[%"PRId64"].filename = %s should equal locations filename = %s\n",
                i, new_locations[i].filename, locations[i].filename);

        
        assert(new_locations[i].forestid == locations[i].forestid);
        assert(new_locations[i].tree_root == locations[i].tree_root);
        assert(new_locations[i].fileid == locations[i].fileid);
        assert(new_locations[i].offset == locations[i].offset);
        assert(new_locations[i].bytes == locations[i].bytes);
        /* fprintf(stderr,"locations[%"PRId64"].bytes = %"PRId64"\n", */
        /*         i,locations[i].bytes); */
    }
    
    /* Check that the preceeding bytes computation is correct */
    {
        int64_t *total_tree_bytes = my_calloc(sizeof(*total_tree_bytes), nfiles);
        for(int64_t i=0;i<ntrees;i++) {
            /* add the number of bytes for tree in each file */
            total_tree_bytes[locations[i].fileid] += locations[i].bytes;
        }
        
        for(int i=0;i<nfiles;i++) {
            XASSERT(total_tree_bytes[i] < inp_file_sizes[i],
                    "Bytes in tree = %"PRId64" must be smaller than file size = %"PRId64"\n",
                    total_tree_bytes[i], inp_file_sizes[i]);
        }
        free(total_tree_bytes);
    }

    
    /* Now assign all trees in the same forest to the same file
       The new fileids goes into new_locations (which is otherwise a copy of locations)
     */
    assign_trees_in_forest_to_same_file(ntrees, locations, new_locations, nfiles, BOX_DIVISIONS);

    /* Now write out both the old and the new struct locations */
    my_snprintf(buffer, MAXLEN, "%s/forests_and_locations_old.list",output_dir);
    write_forests_and_locations(buffer, ntrees, locations);

    /* write new the forests file */
    my_snprintf(buffer,MAXLEN,"%s/forests.list", output_dir);
    unlink(buffer);
    FILE *fp_forests = my_fopen(buffer,"w");
    fprintf(fp_forests, "#TreeRootID ForestID\n");
    for(int64_t i=0;i<ntrees;i++) {
        fprintf(fp_forests, "%"PRId64" %"PRId64"\n",
                locations[i].tree_root, locations[i].forestid);
    }
    fclose(fp_forests);
    
    /* open the locations file*/
    my_snprintf(buffer,MAXLEN,"%s/locations.dat", output_dir);
    unlink(buffer);
    FILE *fp_locations = my_fopen(buffer,"w");
    fprintf(fp_locations, "#TreeRootID FileID Offset Filename\n");

    
    /* copy the headers between the tree_* files */
    /* break when the number of trees is encountered -- should be the first one line that doesn't have a '#' character at front */
    int64_t *tree_header_offsets = my_malloc(sizeof(*tree_header_offsets), nfiles);
    for(int i=0;i<nfiles;i++) {
        /* All of the file pointers have been moved around to figure out the bytes
           -> reposition them at the beginning of the tree_*.dat file
         */
        rewind(tree_inputs[i]);
        
        while(fgets(buffer, MAXLEN, tree_inputs[i]) != NULL) {
            if(buffer[0] != '#') {
                tree_header_offsets[i] = ftello(tree_outputs[i]);
                /* write a place holder for the number of trees in the file.
                   There are 18 X's in the following line, DO NOT CHANGE. 
                 */
                fprintf(tree_outputs[i], "XXXXXXXXXXXXXXXXXX\n"); //For the number of trees                
                break;
            } else {
                fprintf(tree_outputs[i], "%s", buffer);
            }
        }
        tree_outputs_fd_offset[i] = ftello(tree_outputs[i]);
    }

    /* Figure out the offsets and write out a binary file containing the new locations info */
    for(int64_t i=0;i<ntrees;i++) {
        const int tree_bytes_line_size = my_snprintf(buffer, MAXLEN, "#tree %"PRId64"\n", locations[i].tree_root);
        const int64_t bytes_to_write   = locations[i].bytes;
        const int64_t out_fileid = new_locations[i].fileid;
        XASSERT(out_fileid < nfiles,
                "Output fileid = %"PRId64" must be smaller than total number of files = %d\n" ,
                out_fileid, nfiles);
        /* XASSERT(new_locations[i].bytes == bytes_to_write, */
        /*         "new locations bytes = %"PRId64"should be identical to old locations bytes = %"PRId64"\n", */
        /*         new_locations[i].bytes,bytes_to_write); */
        new_locations[i].offset = tree_outputs_fd_offset[out_fileid] + tree_bytes_line_size;
        tree_outputs_fd_offset[out_fileid] += (bytes_to_write + tree_bytes_line_size);
    }


    /* Valgrind complains there is use of uninitialized bytes -> so ditching this binary file output for now */
/* /\* Output the binary locations struct so I can skip over recalculating the bytes *\/ */
    /* { */
    /*     my_snprintf(buffer, MAXLEN, "%s/new_locations.binary",output_dir); */
    /*     FILE *fp = my_fopen(buffer, "w"); */
    /*     /\* fprintf(stderr,"ntrees = %"PRId64"\n",ntrees); *\/ */
    /*     my_fwrite(&ntrees, sizeof(int64_t), 1, fp); */
    /*     const size_t size_of_struct = sizeof(struct locations); */
    /*     /\* fprintf(stderr,"struct size = %zu\n", size_of_struct); *\/ */
    /*     my_fwrite(&size_of_struct, sizeof(size_t), 1, fp); */
    /*     /\* my_fwrite(new_locations, size_of_struct, ntrees, fp); *\/ */
    /*     fclose(fp); */
    /* } */

    /* Write out the combined forests and locations file */
    my_snprintf(buffer, MAXLEN, "%s/forests_and_locations_new.list",output_dir);
    write_forests_and_locations(buffer, ntrees, new_locations);
    
    fprintf(stderr, ANSI_COLOR_MAGENTA"Writing out trees in contiguous order...."ANSI_COLOR_RESET"\n");
    interrupted=0;
    init_my_progressbar(ntrees, &interrupted);

    /* Now copy each one of the trees */
    for(int64_t i=0;i<ntrees;i++) {
        my_progressbar(i, &interrupted);

        const int64_t fileid = locations[i].fileid;
        XASSERT(locations[i].tree_root == new_locations[i].tree_root,
                "locations->tree_root = %"PRId64" must equal new_locations->tree_root = %"PRId64"\n",
                locations[i].tree_root, new_locations[i].tree_root);

        XASSERT(locations[i].forestid == new_locations[i].forestid,
                "locations->forestid = %"PRId64" must equal new_locations->forestid = %"PRId64"\n",
                locations[i].forestid, new_locations[i].forestid);
        
        /* Make sure all output is done using new_locations[i].fileid */
        const int64_t out_fileid = new_locations[i].fileid;
        FILE *out_fp = tree_outputs[out_fileid];
        /* const int tree_bytes_line_size = fprintf(out_fp, "#tree %"PRId64"\n", locations[i].tree_root); */
        fprintf(out_fp, "#tree %"PRId64"\n", locations[i].tree_root);
        fflush(out_fp);

        const int64_t offset          = locations[i].offset;
        const int64_t bytes_to_write  = locations[i].bytes;

        if(bytes_to_write == 0) {
            fprintf(stderr, "Strange! bytes for tree data = %zu should not be 0\n", bytes_to_write);
            continue;
        }
#ifdef USE_FGETS //USE_FGETS -> stdio.h family
#warning using fgets (slowest)
        FILE *in_fp  = tree_inputs[fileid];
        my_fseek(in_fp, (long) offset, SEEK_SET);
        const long actual_offset = ftello(out_fp);
        XASSERT(actual_offset == new_locations[i].offset,
                "actual offset = %ld should equal calculated offset = %"PRId64"\n",
                actual_offset, new_locations[i].offset);
        /* new_locations[i].offset = ftello(out_fp);                 */
        const int64_t bytes_written = copy_bytes_between_two_files(bytes_to_write, in_fp, out_fp);
        
#else //use pread/write
#warning using pread
        int in_fd  = tree_inputs_fd[fileid];
        int out_fd = tree_outputs_fd[out_fileid];
        off_t in_offset = offset;
        const int64_t bytes_written = copy_bytes_with_pread(bytes_to_write, in_fd, out_fd, in_offset);

        /* I have already figured out the offsets */
        /* new_locations[i].offset = tree_outputs_fd_offset[out_fileid] + tree_bytes_line_size; */
#endif//USE_FGETS -> stdio.h family
        
        
        XASSERT(bytes_written == bytes_to_write,
                "bytes_to_write = %zu does not equal bytes_written = %zu\n",
                bytes_to_write, bytes_written);

        /* Update the number of trees in that file */
        tree_counts[out_fileid]++;

        /* write the locations info*/
        const int ii = out_fileid/(BOX_DIVISIONS*BOX_DIVISIONS);
        const int jj = (out_fileid%((int64_t)(BOX_DIVISIONS*BOX_DIVISIONS)))/BOX_DIVISIONS;
        const int kk = out_fileid%((int64_t)BOX_DIVISIONS);
        fprintf(fp_locations, "%"PRId64" %"PRId64" %"PRId64" tree_%d_%d_%d.dat\n",
                new_locations[i].tree_root, out_fileid, new_locations[i].offset, ii, jj, kk);

        /* This line is only required if offsets have not been computed earlier */
        /* tree_outputs_fd_offset[out_fileid] += (bytes_written + tree_bytes_line_size) ; */
    }

    /* fill in the number of trees written per file. the number in the format
     *MUST EXACTLY* match the number of XXX's in the previous place-holder. 
     */
    for(int i=0;i<nfiles;i++) {
        FILE *out_fp = tree_outputs[i];
        fseek(out_fp, tree_header_offsets[i], SEEK_SET);
        fprintf(out_fp, "%-18"PRId64"\n", tree_counts[i]);
    }
    finish_myprogressbar(&interrupted);
    fprintf(stderr, ANSI_COLOR_GREEN "Writing out trees in contiguous order.....done"ANSI_COLOR_RESET"\n\n");


    /* close open file pointers + free memory for file pointers */
    fclose(fp_locations);
    for(int i=0;i<nfiles;i++) {
        fclose(tree_inputs[i]);
        fclose(tree_outputs[i]);
    }
    free(tree_inputs);free(tree_outputs);
    free(tree_inputs_fd);free(tree_outputs_fd);

    /* free other heap allocations */
    free(tree_header_offsets);
    free(tree_outputs_fd_offset);
    free(tree_counts);
    free(inp_file_sizes);
    free(locations);
    free(new_locations);

    gettimeofday(&tend, NULL);
    fprintf(stderr,"Wrote out %"PRId64" trees in contiguous order. Time taken = %0.2g seconds\n",
            ntrees, ADD_DIFF_TIME(tstart, tend));
    
    return EXIT_SUCCESS;
}
    
    
