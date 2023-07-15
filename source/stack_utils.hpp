#ifndef STACK_UTILS_HPP
#define STACK_UTILS_HPP

#include <ap_int.h>
#include <hls_stream.h>

template <typename DATA_T,
          size_t DATA_LENGTH,
          size_t MAX_SOL,
          size_t STREAM_DEPTH,
          size_t STACK_DEPTH>
void stack(
    hls::stream<DATA_T> &in_stream,
    hls::stream<DATA_T> &out_stream,
    hls::stream<bool> &stop_req_stream,
    hls::stream<bool> &overflow_stream)
{
    const DATA_T MASK_NEW_SOLUTION = ~(1UL << (DATA_LENGTH - 1));
    const DATA_T MASK_END_EXTENSION = ~(1UL << (DATA_LENGTH - 2));
    const size_t BIT_NEW_SOLUTION = DATA_LENGTH - 1;
    const size_t BIT_END_EXTENSION = DATA_LENGTH - 2;
    enum STATE_T
    {
        evaluate,
        push,
        pop
    };

    unsigned int top = 0;
    unsigned int max_top = 0;
    STATE_T state{evaluate};
    DATA_T stack[STACK_DEPTH];
    DATA_T curr_sol[MAX_SOL];
    DATA_T curr_sol_length = 0;

    while (true)
    {
        switch (state)
        {

        case evaluate:
        {

            if (in_stream.size() > 0)
            {
                // push solution
                state = push;
            }
            else if ((STREAM_DEPTH - out_stream.size()) > 2 && top > 0)
            {
                // pop solution
                state = pop;
            }
            else
            {
                state = evaluate;
            }

            break;
        }

        case push:
        {
            DATA_T node;
            DATA_T sol_length = in_stream.read();
            // std::cout << "pushing sol of length: " << sol_length << std::endl;

            // Push old solution in reverse order, if any
        PUSH_OLD_SOL_TO_STACK_LOOP:
            for (int g = curr_sol_length - 1; g >= 0; g--)
            {
                stack[top++] = curr_sol[g];
            }
            if (curr_sol_length != 0)
            {
                stack[top++] = curr_sol_length | ~MASK_NEW_SOLUTION;
            }

            // Save new top solution
            curr_sol_length = sol_length & MASK_NEW_SOLUTION;
        PUSH_NEW_TOP_SOL_LOOP:
            for (int g = 0; g < curr_sol_length; g++)
            {
                curr_sol[g] = in_stream.read();
            }

            // Push expansions
            node = in_stream.read();
            stack[top++] = node | ~MASK_END_EXTENSION;
        
        PUSH_NEW_EXTENSIONS:
            while (!node.test(BIT_END_EXTENSION))
            {
                node = in_stream.read();
                stack[top++] = node & MASK_END_EXTENSION;
                if (top == STACK_DEPTH)
                {
                    overflow_stream.write(true);
                    break;
                }
            }

            if (top > max_top){
                max_top = top;
                std::cout << max_top << std::endl;
            }

            state = evaluate;
            break;
        }

        case pop:
        {

            DATA_T node;

            // Write top solution
            out_stream.write(curr_sol_length | ~MASK_NEW_SOLUTION);
        POP_WRITE_OUT_TOP_SOL_LOOP:
            for (int g = 0; g < curr_sol_length; g++)
            {
                out_stream.write(curr_sol[g]);
            }

            // Write expansion
            node = stack[--top];
            out_stream.write(node & MASK_END_EXTENSION);

            // Load new top solution from stack
            if ((node.test(BIT_END_EXTENSION)) && top > 0)
            {
                DATA_T sol_length = stack[--top];
                curr_sol_length = sol_length & MASK_NEW_SOLUTION;
            POP_SOL_FROM_STACK_LOOP:
                for (int g = 0; g < curr_sol_length; g++)
                {
                    curr_sol[g] = stack[--top];
                }
            }

            state = evaluate;
            break;
        }
        }

        bool req_s;
        if (stop_req_stream.read_nb(req_s))
        {
            std::cout << max_top << std::endl;
            break;
        }
    }
}
#endif /*STACK_UTILS_HPP*/