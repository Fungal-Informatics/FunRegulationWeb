import { uuid } from "../helpers/utils";
import { Board, BoardType, SimpleUser } from "../sdkgen/api-generated";
import { QUERY } from "../sql";
import {
	ISelectBoardsSharedWithUserSqlResult,
	ISelectUserOwnedBoardsSqlResult,
} from "../sql/elaborated.query";

async function createBoard(userId: string, name: string, type: BoardType) {
	const newId = uuid();

	if (!name) throw new Error("Nome inválido");

	await QUERY.boards.createBoard(userId, newId, name, type);
}

async function updateBoard(
	userId: string,
	boardId: string,
	name: string,
) {
	const board = await QUERY.boards.getBoard(boardId);

	if (userId !== board?.ownerId) throw new Error("Sem acesso a essa board");

	await QUERY.boards.updateBoard(boardId, name);
}

async function archiveBoard(userId: string, boardId: string) {
	const board = await QUERY.boards.getBoard(boardId);

	if (userId !== board?.ownerId) throw new Error("Sem acesso a essa board");
	if (board?.isArchived) throw new Error("A board já está arquivada");

	await QUERY.boards.archiveBoard(boardId);
}

function transformDbBoard(
	board: ISelectBoardsSharedWithUserSqlResult | ISelectUserOwnedBoardsSqlResult,
	userOwned: boolean,
): Board {
	const types: Record<string, BoardType> = {
		bookmarks: "bookmarks",
		notes: "notes",
		task: "task",
	};

	const shared = (board.sharedWith as unknown) as SimpleUser[];

	return {
		ownedByMe: userOwned,
		sharedWith: shared ?? null,
		type: board.type === null ? null : types[board.type],
		id: board.id,
		name: board.name,
		createdAt: new Date(board.createdAt),
		updatedAt: new Date(board.updatedAt),
		isArchived: board.isArchived,
	};
}

async function getUserBoards(userId: string): Promise<Board[]> {
	// USER OWNED
	const userOwnedBoards = await QUERY.boards.getUserBoards(userId);
	const userOwnedBoardsTreated = userOwnedBoards.map((b) =>
		transformDbBoard(b, true),
	);

	// SHARED WITH USER
	const boardsSharedWithUser = await QUERY.boards.getBoardsSharedWithUser(
		userId,
	);
	const boardSharedWithUserTreated = boardsSharedWithUser.map((b) =>
		transformDbBoard(b, false),
	);

	return [...userOwnedBoardsTreated, ...boardSharedWithUserTreated];
}

async function getBoard(userId: string, boardId: string): Promise<Board> {
	const board = await QUERY.boards.getBoard(boardId);
	if (!board) throw new Error("Board inexistente");

	const shared = (board.sharedWith as unknown) as SimpleUser[] | null;

	if (board.ownerId === userId) {
		return transformDbBoard(board, true);
	} else if (shared !== null && shared.map((b) => b.id).includes(userId)) {
		return transformDbBoard(board, false);
	} else {
		throw new Error("Sem acesso a essa board");
	}
}

async function shareBoard(
	userId: string,
	boardId: string,
	shareWithEmail: string, // another user email
) {
	const targetUser = await QUERY.users.getUserByEmail(shareWithEmail);
	if (!targetUser) throw new Error("Usuario não cadastrado");

	const board = await QUERY.boards.getBoard(boardId);
	if (!board) throw new Error("Board inexistente");

	if (userId !== board.ownerId) throw new Error("Sem acesso a essa board");
	if (board?.isArchived) throw new Error("A board está arquivada");

	await QUERY.boards.createBoardSharing(uuid(), boardId, targetUser.id);
}

export const boards = {
	createBoard,
	updateBoard,
	shareBoard,
	getUserBoards,
	getBoard,
	archiveBoard,
};
