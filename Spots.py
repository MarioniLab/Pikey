import numpy as np
import itertools


def findLocalMaxima(image, thresh):
    '''
    Given some threhsold, find the local maxima of pixel intensities
    in an image.
    '''
    
    max_pixels = image >= thresh
    array_indices = np.where(max_pixels)

    # get local maxima as tuples of row, col indices
    index_tuples = zip(array_indices[0],
                       array_indices[1])
    local_maxima = [q for q in index_tuples]
    
    return(local_maxima)

def findStackLocalMaxima(im_stack, threshold_stack):
    '''
    Iterate over a set of stacked images, with pre-computed
    thresholds for each image to find local maxima
    '''
    maxima_dict = {}
    for im in range(im_stack.shape[2]):
        maxima = findLocalMaxima(im_stack[:,:,im],
                                 threshold_stack[im])
        maxima_dict[str(im)] = maxima

    return(maxima_dict)

def assignBarCodes(pixel_dict, im_stack, stack_thresholds):
    '''
    Assign barcodes to each pixel in the `pixel_dict` given
    it surpasses some threshold within each hybridisation.
    
    '''

    for pos in pixel_dict.keys():
        ix = pos[0]
        jx = pos[1]
        bin_vec = pixel_dict[pos]
    
        # does this position fall into the maximal set
        # in any image?
        max_pixels = im_stack[ix, jx, :] >= stack_thresholds
        # get images to set to 1 bit
        array_indices = np.where(max_pixels)
        bin_vec[array_indices] = 1
    
        signal_vec = im_stack[ix, jx, :]
        pixel_dict[pos] = tuple(bin_vec)
    
    return(pixel_dict)
    

def findExactMatches(pixel_dict, gene_codes):
    '''
    Find the exact pixel matches in a code book
    '''
    exact_matches = set()
    for pixel in pixel_dict.keys():
        pix_vec = np.asarray(pixel_dict[pixel])
    
        # find the genes with an exact match, i.e. error free
        matched = tuple(set([tuple(pixel_dict[pixel])]).intersection(gene_codes))
        if len(matched) == 1:
            exact_matches.add(pixel)

    return(exact_matches)

def selectSignalThreshold(image_stack, code_book):
    '''
    Select the threshold at which to call local maxima pixels.
    This is taken as the hybridisation-specific pixel intensities
    at which the number of exact barcode matches is maximised.
    
    NB: In which situations would this utterly fail?  i.e. if
    there is no signal whatsoever?
    
    Args
    ----
    image_stack: np.array
        Numpy array with shape (x, y, N_hyb), where `N_hyb` is
        the number of hybridisations
        
    code_book: dict or list
        mapping of gene barcodes to gene IDs, or list of gene
        barcodes (IDs not required)
        
    Returns
    -------
    thresholds: np.array
        Numpy array of pixel intensity thresholds at which to call
        local maximal intensities (i.e. spots)
    '''
    
    # set a random seed?
    #np.random.seed(42)
    np.random.seed()
    n_hyb = image_stack.shape[2]
    spaces = 10
    n_exact = np.zeros(spaces)
    #poss_thresh = np.linspace(0.01, 0.99, num=spaces)
    max_stack = np.max(image_stack.max(axis=0), axis=0)

    #max_array = np.zeros((n_hyb, spaces))

    if type(code_book) == dict:
        gene_codes = sorted(code_book.keys())
    elif type(code_book) == list:
        gene_codes = code_book

    # need to enumerate all possible combinations of thresholds
    # start with a small number and build up
    # iterate over all possible combinations
    # if the new number of exact matches is bigger than the current,
    # store it, else discard and move on. Don't need to store all of them.
    # what about randomly sampling from the thresholds, for X iterations
    # to find the best set?
    # probably not guaranteed to find the global maxiam??
    
    max_iters = 100
    current_max = 0
    current_comb = np.full((n_hyb, ), 1.0)
    
    it = 0
    while it < max_iters:
    #for it in range(max_iters):
        # randomly sample from threhsolds
        # can we reject any out of hand, i.e all > 0.9 ??
        poss_thresh = np.random.random_sample(n_hyb)
        
        # check that > 10% of thresholds are not higher than the current best
        check = poss_thresh > current_comb
        if np.sum(check)/float(n_hyb) > 0.1:
            pass
        else:
            # this proportion is probably different for each hyb
            # need to find %threshold for each hyb that maximises
            # the number of exact matches
            # this will make it significantly slower
        
            thresh_stack = poss_thresh * max_stack
            maxima_dict = findStackLocalMaxima(image_stack, thresh_stack)
        
            # take the union of all maxima
            union_maxima = set([u for u in itertools.chain.from_iterable(maxima_dict.values())])
            n_maxima = len(union_maxima)
        
            # make a dictionary to store the pixels and their binary arrays
            pixel_dict = dict(zip(union_maxima, (np.zeros((n_hyb,)) for i in range(0, n_maxima))))
            
            # assign the barcodes to each pixel given the thresholds on each
            # hybridisation
            assigned_dict = assignBarCodes(pixel_dict, image_stack, thresh_stack)
        
            # find the pixels where there is an exact match to a gene barcode
            perfect_match = findExactMatches(assigned_dict, gene_codes)
        
            # check if this is a new maxima, if not discard
            n_exact = len(perfect_match)
        
            if n_exact > current_max:
                current_max = n_exact
                current_comb = thresh_stack
            else:
                pass
            it += 1

        return(current_max, current_comb)
    
    
def selectSignalThresholdBinary(image_stack, code_book):
    '''
    Select the threshold at which to call local maxima pixels.
    This is taken as the hybridisation-specific pixel intensities
    at which the number of exact barcode matches is maximised.
    
    NB: In which situations would this utterly fail?  i.e. if
    there is no signal whatsoever?
    
    Use a binary search to find the optimal threshold for each
    hybridisation (?)
    
    Args
    ----
    image_stack: np.array
        Numpy array with shape (x, y, N_hyb), where `N_hyb` is
        the number of hybridisations
        
    code_book: dict or list
        mapping of gene barcodes to gene IDs, or list of gene
        barcodes (IDs not required)
        
    Returns
    -------
    thresholds: np.array
        Numpy array of pixel intensity thresholds at which to call
        local maximal intensities (i.e. spots)
    '''
    
    # set a random seed?
    #np.random.seed(42)
    np.random.seed()
    n_hyb = image_stack.shape[2]
    #poss_thresh = np.linspace(0.01, 0.99, num=spaces)
    max_stack = np.max(image_stack.max(axis=0), axis=0)

    #max_array = np.zeros((n_hyb, spaces))

    if type(code_book) == dict:
        gene_codes = sorted(code_book.keys())
    elif type(code_book) == list:
        gene_codes = code_book

    # need to enumerate all possible combinations of thresholds
    # start with a small number and build up
    # iterate over all possible combinations
    # if the new number of exact matches is bigger than the current,
    # store it, else discard and move on. Don't need to store all of them.
    # what about randomly sampling from the thresholds, for X iterations
    # to find the best set?
    # probably not guaranteed to find the global maxiam??
    
    # this is going to be slow!!
    # look into optimising this somewhat
    
    max_iters = 100
    current_comb = np.full((n_hyb, ), 0.5)
    
    # initialize with the 0.5 values
    thresh_stack = current_comb * max_stack
    maxima_dict = findStackLocalMaxima(image_stack, thresh_stack)
        
    # take the union of all maxima
    union_maxima = set([u for u in itertools.chain.from_iterable(maxima_dict.values())])
    n_maxima = len(union_maxima)
        
    # make a dictionary to store the pixels and their binary arrays
    pixel_dict = dict(zip(union_maxima, (np.zeros((n_hyb,)) for i in range(0, n_maxima))))
            
    # assign the barcodes to each pixel given the thresholds on each
    # hybridisation
    assigned_dict = assignBarCodes(pixel_dict, image_stack, thresh_stack)
        
    # find the pixels where there is an exact match to a gene barcode
    perfect_match = findExactMatches(assigned_dict, gene_codes)
        
    # check if this is a new maxima, if not discard
    current_max = len(perfect_match)
    old_max = current_max
    old_comb = current_comb
    
    for it in range(max_iters):
        if current_max > old_max:
            old_max = current_max
            old_comb = current_comb
        elif current_max < old_max:
            break
        else:
            pass

        # do a binary search on each hybridisation
        # to find the threshold that maximises the number of exact matches
        # can we reject any out of hand, i.e all > 0.9 ??
        for hx in range(n_hyb):
            curr_val = current_comb[hx]
            low_val = curr_val/2.0
            low_comb = current_comb.copy()
            low_comb[hx] = low_val
                    
            # this proportion is probably different for each hyb
            # need to find %threshold for each hyb that maximises
            # the number of exact matches
            # this will make it significantly slower
        
            thresh_stack = low_comb * max_stack
            maxima_dict = findStackLocalMaxima(image_stack, thresh_stack)
        
            # take the union of all maxima
            union_maxima = set([u for u in itertools.chain.from_iterable(maxima_dict.values())])
            n_maxima = len(union_maxima)
        
            # make a dictionary to store the pixels and their binary arrays
            pixel_dict = dict(zip(union_maxima, (np.zeros((n_hyb,)) for i in range(0, n_maxima))))
            
            # assign the barcodes to each pixel given the thresholds on each
            # hybridisation
            assigned_dict = assignBarCodes(pixel_dict, image_stack, thresh_stack)
        
            # find the pixels where there is an exact match to a gene barcode
            perfect_match = findExactMatches(assigned_dict, gene_codes)
        
            # check if this is a new maxima, if not discard
            n_exact = len(perfect_match)
        
            if n_exact > current_max:
                current_comb = low_comb
                current_max = n_exact
            else:
                hi_val = (curr_val + low_val)/2.0
                hi_comb = current_comb.copy()
                hi_comb[hx] = hi_val
                
                thresh_stack = hi_comb * max_stack
                maxima_dict = findStackLocalMaxima(image_stack, thresh_stack)
        
                # take the union of all maxima
                union_maxima = set([u for u in itertools.chain.from_iterable(maxima_dict.values())])
                n_maxima = len(union_maxima)
        
                # make a dictionary to store the pixels and their binary arrays
                pixel_dict = dict(zip(union_maxima, (np.zeros((n_hyb,)) for i in range(0, n_maxima))))
            
                # assign the barcodes to each pixel given the thresholds on each
                # hybridisation
                assigned_dict = assignBarCodes(pixel_dict, image_stack, thresh_stack)
        
                # find the pixels where there is an exact match to a gene barcode
                perfect_match = findExactMatches(assigned_dict, gene_codes)
        
                # check if this is a new maxima, if not discard
                n_exact = len(perfect_match)
                
                current_max = n_exact
                current_comb = hi_comb
        
    return(old_max, old_comb)
