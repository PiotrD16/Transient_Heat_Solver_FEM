#include "../header/gaussLegendre.h"

double Integral::calculateIntegralOneDim(double (*func)(double), const int points) const {
    double result = 0.0;

    switch (points) {
        case 1: {
            for (int i = 0; i < schema1N.getElements(); ++i)
                result += schema1N.getWeights()[i] * func(schema1N.getNodes()[i]);
        } break;
        case 2: {
            for (int i = 0; i < schema2N.getElements(); i++)
                result += schema2N.getWeights()[i] * func(schema2N.getNodes()[i]);
        }break;
        case 3: {
            for (int i = 0; i < schema3N.getElements(); i++)
                result += schema3N.getWeights()[i] * func(schema3N.getNodes()[i]);
        }break;
        case 4: {
            for (int i = 0; i < schema4N.getElements(); i++)
                result += schema4N.getWeights()[i] * func(schema4N.getNodes()[i]);
        } break;
        default: {
            std::cout << "You can pass 1, 2, 3 or 4 points!\n";
            return NAN;
        }
    }
    return result;
}

double Integral::calulateIntegralTwoDim(double (*func)(double ,double), int points) const {
    double result = 0.0;
    const GaussLegendreSchema* schema = nullptr;

    switch (points) {
        case 1: schema = &schema1N; break;
        case 2: schema = &schema2N; break;
        case 3: schema = &schema3N; break;
        case 4: schema = &schema4N; break;
        default: std::cout << "You can pass 1, 2, 3 or 4 points!\n"; break;
    }

    if (schema == nullptr) {
        std::cout << "Schema is null, exiting...!\n";
        exit(1);
    }

    for (int i = 0; i < schema->getElements(); i++) {
        for (int j = 0; j < schema->getElements(); j++) {
            result += schema->getWeights()[i] * schema->getWeights()[j] * func(schema->getNodes()[i], schema->getNodes()[j]);
        }
    }

    return result;
}