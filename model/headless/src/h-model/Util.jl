module Util

function get_black_pixel_coordinates_scaled(image_path, target_width, target_height; threshold = 0.5)
    img_org = load(image_path)
    img = Gray.(img_org)

    height, width = size(img)
    black_pixel_coordinates = Set{Tuple{Int, Int}}()

    x_scale = target_width / width
    y_scale = target_height / height

    for y in 1:height
        for x in 1:width
            pixel_value = img[y, x]
            if pixel_value <= Gray(threshold)
                scaled_x = round(Int, x * x_scale)
                scaled_y = round(Int, y * y_scale)
                push!(black_pixel_coordinates, (scaled_x + 1, scaled_y + 1))
            end
        end
    end

    return [CartesianIndex(idx) for idx in black_pixel_coordinates]
end

end # module
