#include "ranges.h"
#include "geodesy.h"

void test_overlapping(const SphericalModel *model);
void test_disjoint(const SphericalModel *model);
void test_range(const SphericalModel *model);

int main() {

    // Build a new model
    SphericalModel *model = loadSphericalModel("sph_med_500.bin", 500);
    
    // And create a range object
    test_range(model);
    test_disjoint(model);
    test_overlapping(model);

    return 0;
}

void test_overlapping(const SphericalModel *model) {

    /**========================================================================
     *!                           Test overlapping ranges
     *========================================================================**/
    printf("=================== Overlapping ranges ==============\n");
    range *r0_3 = SphericalModelExtractRange(model, 0, 3);
    range *r2_4 = SphericalModelExtractRange(model, 2, 4);
    range *r0_4 = SphericalModelExtractRange(model, 0, 4);
    // range *r1_4 = SphericalModelExtractRange(model, 1, 4);

    range *r0_4_add = ranges_union(r0_3, r2_4);

    assert(ranges_equal(r0_4_add, r0_4));

    // Non zero start
    range *r15_18 = SphericalModelExtractRange(model, 15, 18);
    range *r19_20 = SphericalModelExtractRange(model, 19, 20);
    range *r15_20 = SphericalModelExtractRange(model, 15, 20);
    range *r15_20_add = ranges_union(r15_18, r19_20);


    assert(ranges_equal(r15_20_add, r15_20));
    printf("r15_20 and r15_20_add are equal\n");

    range *ra = SphericalModelExtractRange(model, 53, 201);
    range *ra1 = SphericalModelExtractRange(model, 202, 202);
    range *r53_202 = SphericalModelExtractRange(model, 53, 202);
    range *r53_202_add = ranges_union(ra, ra1);

    assert(r53_202->card == r53_202_add->card);


    // Vector_print_head_d(r53_202->C_lm, 20);
    // Vector_print_head_d(r53_202_add->C_lm, 20);
    // // Are the first 100 elements the same at least??
    // // the first 100 elements are legit
    // // for (int i = 0; i < r53_202->card / 2; i++) {
    // for (int i = 0; i < 50 / 2; i++) {
    //     assert(r53_202->C_lm->data[i] == r53_202_add->C_lm->data[i]);
    //     assert(r53_202->S_lm->data[i] == r53_202_add->S_lm->data[i]);
    // }


    // // assert(ranges_equal(r53_202, r53_202_add));

    // range *rb = SphericalModelExtractRange(model, 150, 204);

    // range *r53_204 = SphericalModelExtractRange(model, 53, 204);
    // range *r53_add = ranges_union(ra, rb);

    // assert(ranges_equal(r53_204, r53_204));
    // assert(ranges_equal(r53_add, r53_add));
    // // assert(ranges_equal(r53_add, r53_204));
    // printRange(r53_add);


}
void test_disjoint(const SphericalModel *model) {

    /**========================================================================
     *!                           Test disjoint ranges
     *========================================================================**/
    printf("=================== Disjoint ranges ==============\n");
    range *low = SphericalModelExtractRange(model, 0, 3);
    range *high = SphericalModelExtractRange(model, 5, 6);

    range *low_high = ranges_union(low, high);

    printRange(low_high);

    // Verify that the values for C_4^* and S_4^* are 0
    for (int i = 0; i < 5; i++) {
        assert(low_high->C_lm->data[i + 10] == 0);
        assert(low_high->S_lm->data[i + 10] == 0);
    }


}
void test_range(const SphericalModel *model) {

        range *r = as_range(model);

    assert(Matrix_size_d(r->C_lm) == r->card / 2);
    assert(Matrix_size_d(r->S_lm) == r->card / 2);

    range *sub = SphericalModelExtractRange(model, 0, 1);

    assert(sub->card == 6);

    // Assert that the elements stored in sub are the same as the first six elements in r

    for (int i = 0; i < 3; i++) {
        assert(sub->C_lm->data[i] == r->C_lm->data[i]);
        assert(sub->S_lm->data[i] == r->S_lm->data[i]);
    }

    printRange(r);

    Matrix_print_d(sub->C_lm);
    Matrix_print_d(sub->S_lm);


    range *r0_4 = SphericalModelExtractRange(model, 0, 4);
    range *r1_4 = SphericalModelExtractRange(model, 1, 4);
    range *r0_0 = SphericalModelExtractRange(model, 0, 0);

    // Matrix_print_d(r0_0->C_lm);
    // Matrix_print_d(r0_0->S_lm);

    range *r0_4_add = ranges_union(r0_0, r1_4);

    printRange(r0_0);
    printf(" + ");
    printRange(r1_4);
    printf(" = ");
    printRange(r0_4_add);

    Matrix_print_d(r0_4_add->C_lm);
    Matrix_print_d(r0_4->C_lm);

    Matrix_print_d(r0_0->S_lm);
    Matrix_print_d(r1_4->S_lm);

    Matrix_print_d(r0_4_add->S_lm);
    Matrix_print_d(r0_4->S_lm);

    assert(ranges_equal(r0_4, r0_4_add));
}