from abc import ABC
from numpy.random import default_rng

class Noise(ABC):
    def sample(self, rng=default_rng()):
        raise NotImplementedError

    def transform(self, input):
        raise NotImplementedError

    def chain(self, other):
        return ChainedNoise(self, other)
        

class ChainedNoise(Noise):
    def __init__(self, sampler, transformer):
        self.sampler = sampler,
        self.transformer = transformer

    def sample(self, rng=default_rng()):
        error = self.sampler.sample(rng)
        return self.transformer.transform(error)


def chain_noises(*models):
    model = models[0]
    for m in models[1:]:
        model = model.chain(m)
    return model
