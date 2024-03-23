

function custom_gradient_descent(Γ)
    v = flatten(Γ)
    p = 1
    grad = Calculus.gradient(v_action, v)
    grad = (flatten ∘ project ∘ emboss)(grad)

    nm = norm(grad)

    old = v_action(v)
    step = p
    iter = 0

    while true
        v_new = (flatten ∘ project ∘ emboss)(v - step * grad)

        new = v_action(v_new)
        if (new > old)
            step /= 2
            continue
        end

        grad_new = (flatten ∘ project ∘ emboss)(Calculus.gradient(v_action, v_new))

        if (iter % 5 == 0)
            cb(v_new)
        end
        iter += 1
        old, v, grad = new, v_new, grad_new

        nm = norm(grad)
        println(step, " ", new, "      ", nm)

        step = min(p, step * 2)
        iter += 1
        if nm < 1e-6 || iter > 1e6 || step < 1e-1 break end
    end
end
