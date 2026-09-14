function [best, ranking] = simplifyCompact(expr, steps)
%SIMPLIFYCOMPACT Busca una representación compacta de una expresión simbólica.
%
%   best = simplifyCompact(expr)
%   best = simplifyCompact(expr, steps)
%   [best, ranking] = simplifyCompact(...)
%
%   La función genera distintas formas algebraicamente equivalentes
%   y selecciona la representación textual más corta.
%
%   Requiere Symbolic Math Toolbox.

    if nargin < 2
        steps = 100;
    end

    if ~isa(expr, 'sym')
        error('La entrada debe ser una expresión simbólica.');
    end

    if ~isscalar(expr)
        error('La versión actual está diseñada para expresiones escalares.');
    end

    candidates = sym.empty(0,1);
    methods    = strings(0,1);

    % -------------------------------------------------------------
    % Expresión original
    % -------------------------------------------------------------

    addCandidate(expr, "Original");


    % -------------------------------------------------------------
    % simplify(): devuelve formas equivalentes
    % -------------------------------------------------------------

    try
        s = simplify(expr, Steps=steps, All=true);

        for k = 1:numel(s)
            addCandidate(s(k), "simplify");
        end
    catch
    end


    % -------------------------------------------------------------
    % Fracción racional
    % -------------------------------------------------------------

    try
        s = simplifyFraction(expr);
        s = simplify(s, Steps=steps);
        addCandidate(s, "simplifyFraction");
    catch
    end


    % -------------------------------------------------------------
    % Expandir antes de simplificar
    % -------------------------------------------------------------

    try
        s = expand(expr);
        s = simplify(s, Steps=steps);
        addCandidate(s, "expand -> simplify");
    catch
    end


    % -------------------------------------------------------------
    % Combinar estructuras algebraicas
    % -------------------------------------------------------------

    try
        s = combine(expr);
        s = simplify(s, Steps=steps);
        addCandidate(s, "combine -> simplify");
    catch
    end


    % -------------------------------------------------------------
    % Combinación trigonométrica
    % -------------------------------------------------------------

    try
        s = combine(expr, 'sincos');
        s = simplify(s, Steps=steps);
        addCandidate(s, "combine(sincos) -> simplify");
    catch
    end


    % -------------------------------------------------------------
    % Forma factorizada
    % factor() devuelve los factores por separado.
    % Se reconstruye el producto con prod().
    % -------------------------------------------------------------

    try
        factors = factor(expr);
        s = prod(factors);
        s = simplify(s, Steps=steps);
        addCandidate(s, "factor -> simplify");
    catch
    end


    % -------------------------------------------------------------
    % Collect respecto de cada variable
    % -------------------------------------------------------------

    vars = symvar(expr);

    for k = 1:numel(vars)

        try
            s = collect(expr, vars(k));
            s = simplify(s, Steps=steps);

            addCandidate( ...
                s, ...
                "collect(" + string(vars(k)) + ") -> simplify" ...
            );
        catch
        end

    end


    % -------------------------------------------------------------
    % Eliminar expresiones duplicadas
    % -------------------------------------------------------------

    repr = string(candidates);

    [~, idx] = unique(repr, 'stable');

    candidates = candidates(idx);
    methods    = methods(idx);
    repr       = repr(idx);


    % -------------------------------------------------------------
    % Métrica de complejidad
    %
    % Aquí se usa simplemente la cantidad de caracteres.
    % Cuanto menor sea, más compacta se considera la expresión.
    % -------------------------------------------------------------

    score = strlength(repr);

    [score, order] = sort(score);

    candidates = candidates(order);
    methods    = methods(order);
    repr       = repr(order);


    % -------------------------------------------------------------
    % Mejor resultado
    % -------------------------------------------------------------

    best = candidates(1);


    % -------------------------------------------------------------
    % Tabla de resultados
    % -------------------------------------------------------------

    ranking = table( ...
        methods, ...
        repr, ...
        score, ...
        'VariableNames', {'Metodo', 'Expresion', 'Complejidad'} ...
    );


    % =============================================================
    % Función auxiliar
    % =============================================================

    function addCandidate(candidate, method)

        if isscalar(candidate)
            candidates(end+1,1) = candidate;
            methods(end+1,1)    = method;
        end

    end

end