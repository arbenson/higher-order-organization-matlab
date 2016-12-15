function z = NMI(classes, clusters)
%NMI - calculates normalized mutual information to evaluate clustering
% z = NMI(classes, clusters) where classes are the true classes
% and clusters is the clustering assignment.

assert(numel(classes) == numel(clusters));

uw = unique(classes);
uc = unique(clusters);

Cw = counts(classes, uw);
Cc = counts(clusters, uc);

I = 0;
N = numel(classes);

for k = 1:length(uw)
    for j = 1:length(uc)
        common = intersect(find(classes == uw(k)), find(clusters == uc(j)));
        nc = numel(common);
        if nc == 0
            continue;
        end
        I = I + (nc / N) * log(N * nc / Cw(k) / Cc(j));
    end
end

z = 2 * I / (entropy(classes) + entropy(clusters));
end

function C = counts(x, ux)
    N = numel(ux);
    C = zeros(N, 1);
    for i = 1:length(ux)
        C(i) = sum(x == ux(i));
    end
end

function H = entropy(x)
    ux = unique(x);
    H = 0;
    N = numel(x);
    for i = 1:length(ux)
        nuxi = sum(x == ux(i));
        H = H - (nuxi / N) * log(nuxi / N);
    end
end
