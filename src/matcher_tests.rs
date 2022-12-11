#[cfg(test)]
mod tests {
    use crate::{
        border_matcher::match_borders,
        parsed_puzzles::ParsedPuzzles,
        utils::{load_image_from_path, Side},
    };

    #[test]
    pub fn best_matches() {
        println!("Hello?");
        let color_image = load_image_from_path("img/crop.jpg").unwrap();
        let parsed_puzzles = ParsedPuzzles::new(&color_image);

        let test = |figure_id: usize, b_id: usize, best_figure: usize| {
            let mut options = vec![];

            let figure = &parsed_puzzles.figures[figure_id];
            for other_figure_id in 0..parsed_puzzles.figures.len() {
                if other_figure_id == figure_id {
                    continue;
                }
                let other_figure = &parsed_puzzles.figures[other_figure_id];
                if other_figure.is_good_puzzle() {
                    for j in 0..4 {
                        if let Some(result) = match_borders(
                            &parsed_puzzles,
                            Side {
                                fig: figure_id,
                                side: b_id,
                            },
                            Side {
                                fig: other_figure_id,
                                side: j,
                            },
                        ) {
                            options.push(result);
                        }
                    }
                }
            }
            options.sort_by(|a, b| a.score.total_cmp(&b.score));
            let position = options
                .iter()
                .position(|o| o.rhs_id == best_figure)
                .unwrap();
            eprintln!(
                "{}my_figure = {figure_id}, rot = {b_id}, need = {best_figure}, place = {position}",
                if position == 0 { "" } else { "!!! " }
            );
        };
        eprintln!("Parsed! Start tests.");
        test(39, 0, 9);
        test(39, 3, 16);
        test(278, 0, 320); // absolutely sure
        test(278, 2, 243); // absolutely sure
    }
}
